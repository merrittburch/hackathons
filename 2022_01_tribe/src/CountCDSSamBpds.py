@file:DependsOn("com.github.samtools:htsjdk:2.19.0")
@file:DependsOn("net.maizegenetics:tassel:5.2.60")

import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.ValidationStringency
import net.maizegenetics.util.Utils
import java.io.File
import net.maizegenetics.util.Sizeof

data class TranscriptBpAlignmentStats(val xOrEqArray: MutableList<Int>, val eqArray:MutableList<Int>, var nCounts:Int = 0, var dCounts:Int = 0, var numAlignments:Int = 0 )

enum class ExonOperator {
    L, C, T
}
val exonOperatorMap = mapOf<Char,ExonOperator>('L' to ExonOperator.L,'C' to ExonOperator.C,'T' to ExonOperator.T)
data class ExonBoundary(val length: Int, val operator: ExonOperator)


/**
 * Function to parse out the exon boundaries from the transcript name.
 * They take the form: transcript_1234:100L24C-12C-90C10T
 */
fun parseExonBoundaries(transcriptName: String) : List<List<ExonBoundary>> {
    val exonStrings = transcriptName.split(":")[1].split("-")

    return exonStrings.map { convertExonStringToBoundary(it) }

}

/**
 * Function to take a single exon and parse out the operator and the length of that operation.
 * It needs to have a running sub list of the number characters until it sees an operator.
 * Then it will convert the temp String into an Int and save out the boundary.
 */
fun convertExonStringToBoundary(exonString: String) : List<ExonBoundary> {
    var runningCount = StringBuilder()
    val exonBoundaries = mutableListOf<ExonBoundary>()
    for(element in exonString) {
        //Check to see if we hit an operator
        if(element in exonOperatorMap.keys) {
            //if so, turn the runningCount into an Int and add the boundary to the list.
            val countString = runningCount.toString()
            check(!countString.isEmpty()) {"Error parsing exon name.  No counts for operator"}

            exonBoundaries.add(ExonBoundary(countString.toInt(), exonOperatorMap[element]!!))
            runningCount.clear()
        }
        else {
            runningCount.append("$element")
        }
    }

    return exonBoundaries
}

/**
 * Function from original method to count each cds bp and aggregate across all the sam files
 */
fun getBpSummary(array: MutableList<Int>): List<Int> {
    val countList = mutableListOf<Int>(0,0,0)

    for(i in array.indices) {
        countList[i%3] += array[i]
    }
    return countList
}

/**
 * Function to parse the SAM records and from the CIGAR string determine mapping and eq bps
 * This will make an array for both the Mapping Bps and the EQ bps.
 */
fun computeSamBpArrays(samFile: File, minAlignmentPercent : Double = 0.9, cdsConservPercentage : MutableMap<String,Map<String,Double>> = mutableMapOf()) : Map<String,TranscriptBpAlignmentStats> {
    println("Processing: ${samFile.name}")
    val reader = SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)
        .open(samFile)

    val samIterator = reader.iterator()

    val statMap = mutableMapOf<String, TranscriptBpAlignmentStats>()
    val numAlignmentCountMap = mutableMapOf<String,Int>()
    val cdsConserveMapByTranscript = mutableMapOf<String,Double>()
    while(samIterator.hasNext()) {
        val currentRecord = samIterator.next()
        //Skip over if its unmapped
        if(currentRecord.readUnmappedFlag) {
            continue
        }

//Uncomment for su1 test
//        if(currentRecord.readName != "Zm00001eb174590_T001:224L363C-186C-93C-153C-77C-142C-108C-86C-72C-127C-87C-93C-159C-84C-81C-83C-118C-258C288T") {
//            continue
//        }

        //Increment the Number of alignments for a given Transcript
        numAlignmentCountMap[currentRecord.readName] = (numAlignmentCountMap[currentRecord.readName]?:0)+1

        //If its not secondary, we can build the alignmentBp arrays
        if(!statMap.containsKey(currentRecord.readName) && !currentRecord.isSecondaryOrSupplementary) {
            //Check min
            val exonBoundaries = parseExonBoundaries(currentRecord.readName)
            val cdsBoundaries = computeCDSPositions(exonBoundaries)
            val stats = buildTranscriptBpAlignmentStats(currentRecord)

//            println(currentRecord.readName)
//            println(exonBoundaries)
//            println(stats.xOrEqArray.size)
            val numMapping = stats.xOrEqArray.slice(cdsBoundaries.first .. cdsBoundaries.second).sum()

            val alignmentPercentage = (numMapping.toDouble()/(cdsBoundaries.second - cdsBoundaries.first + 1))
            cdsConserveMapByTranscript[currentRecord.readName] = alignmentPercentage
            if(cdsBoundaries.second != cdsBoundaries.first &&  alignmentPercentage >= minAlignmentPercent ) {
                statMap[currentRecord.readName] = stats
            }

        }

    }

    cdsConservPercentage[samFile.name] = cdsConserveMapByTranscript.toMap()

    //Add alignment counters to statmaps
    for(key in statMap.keys) {
        statMap[key]!!.numAlignments = numAlignmentCountMap[key]?:0
    }
    return statMap
}

fun computeCDSPositions(exonBoundaries:List<List<ExonBoundary>>) : Pair<Int,Int> {
    val sizesOfOperators = exonBoundaries.flatten()
        .groupBy { it.operator }
        .map { Pair(it.key, it.value.map { currValue -> currValue.length }.sum()) }
        .toMap()

    val leaderSize = sizesOfOperators[ExonOperator.L]?:0

    val cdsSize = sizesOfOperators[ExonOperator.C]?:0

    return Pair(leaderSize, leaderSize + cdsSize - 1)
}

/**
 * Function that actually creates the xOrEQ array and the eqArray.  It also counts the number of N and D bps
 */
fun buildTranscriptBpAlignmentStats(samRecord: SAMRecord) : TranscriptBpAlignmentStats {
    val xOrEqArray = Array<Int>(samRecord.readLength) { 0 }
    val eqArray = Array<Int>(samRecord.readLength) { 0 }

    var nCounts = 0
    var dCounts = 0

    var currentBp = 0
    val cigarElements = samRecord.cigar.cigarElements

    //Loop through the cigarElements
    for(element in cigarElements) {
        val operator = element.operator
        val count = element.length

        if(operator== CigarOperator.N) {
            nCounts+=count
        }
        else if(operator==CigarOperator.D) {
            dCounts+=count
        }
        //Check to see if consumes query
        else if(operator.consumesReadBases()) {
            //If it consumes read bases, we can walk through the length of the CIGAR operator and set the position as a 1
            for(index in 0 until count) {
                if(operator.isAlignment) {
                    xOrEqArray[currentBp] = 1
                }

                if (operator==CigarOperator.EQ) {
                    eqArray[currentBp] = 1
                }

                currentBp++
            }
        }

    }

    //Check to see if it was reversed during alignment.  If so we need to flip our arrays.
    return if(samRecord.readNegativeStrandFlag) {
        TranscriptBpAlignmentStats(xOrEqArray.reversed().toMutableList(), eqArray.reversed().toMutableList(), nCounts, dCounts,1)
    }
    else {
        TranscriptBpAlignmentStats(xOrEqArray.toMutableList(), eqArray.toMutableList(),nCounts, dCounts,1)
    }

}

/**
 * Funtion to aggregate multiple SAM file's TranscriptBpAlignmentStats easily.
 * Basically it will check to make sure the lengths of the two arrays are the right size and then will sum the values together.
 */
fun addToAggregate(aggregateMap:MutableMap<String, TranscriptBpAlignmentStats>, currentMap:Map<String, TranscriptBpAlignmentStats>) {
    for((currentTranscript, alignmentStats) in currentMap) {
        if(aggregateMap.containsKey(currentTranscript)) {
            check(aggregateMap[currentTranscript]!!.xOrEqArray.size == alignmentStats.xOrEqArray.size) {"XorEQ arrays are different sizes.  ${currentTranscript}"}
            check(aggregateMap[currentTranscript]!!.eqArray.size == alignmentStats.eqArray.size) {"EQ arrays are different sizes.  ${currentTranscript}"}

            for(i in alignmentStats.xOrEqArray.indices) {
                aggregateMap[currentTranscript]!!.xOrEqArray[i] += alignmentStats.xOrEqArray[i]
                aggregateMap[currentTranscript]!!.eqArray[i] += alignmentStats.eqArray[i]
            }

            aggregateMap[currentTranscript]!!.nCounts += alignmentStats.nCounts
            aggregateMap[currentTranscript]!!.dCounts += alignmentStats.dCounts
            aggregateMap[currentTranscript]!!.numAlignments += alignmentStats.numAlignments
        }
        else {
            aggregateMap[currentTranscript] = alignmentStats
        }
    }
}


fun countBpsInCDSSamDir(samDir: String, outputFile: String, useOriginalMethod:Boolean,useAllSams:Boolean = false, minCDSPercentage : Double = .90, outputCDSPercentageFile: String) {
    val aggregateMap = mutableMapOf<String, TranscriptBpAlignmentStats>()
    val samFiles = File(samDir).walk()
        .filter { it.isFile }
        .filter { !it.isHidden }

    val numSams = if(useAllSams) samFiles.toList().count() else -1

    val cdsCounter = mutableMapOf<String,Map<String,Double>>()
    val samCountByTranscript = samFiles
        .map { computeSamBpArrays(it,minCDSPercentage, cdsCounter) }
        .map {
            addToAggregate(aggregateMap, it)
            it.keys
        }
        .flatten()
        .groupingBy { it }
        .eachCount()



    Utils.getBufferedWriter(outputCDSPercentageFile).use { output ->
        val samFiles = cdsCounter.keys.sorted()
        output.write("TranscriptName\t${samFiles.joinToString("\t")}\n")

        for (transcript in samCountByTranscript.keys.sorted()) {
            output.write("${transcript}\t${samFiles.map { cdsCounter[it]?.get(transcript)?:0.0 }.joinToString("\t")}\n")
        }
    }



    //Total % bp1+bp2 and bp3 for each transcript
    //transcript\tlenTranscript\tnumSams\tbp1CountM\tbp2CountM\tbp3CountM\tbp1CountEQ\tbp2CountEQ\tbp3CountEQ

    if(useOriginalMethod) {
        outputCountsOriginalMethod(outputFile, samCountByTranscript, aggregateMap)
    }
    else {
        outputCountsFullTranscript(outputFile, samCountByTranscript,numSams, aggregateMap)
    }

}

/**
 * Original method to count the bps for CDS based alignments.
 */
fun outputCountsOriginalMethod(
    outputFile: String,
    samCountByTranscript: Map<String, Int>,
    aggregateMap: MutableMap<String, TranscriptBpAlignmentStats>
) {
    Utils.getBufferedWriter(outputFile).use { output ->
        output.write(
            "TranscriptName\tTranscriptLength\tNumSams\t" +
                    "bp1CountM\tbp2CountM\tbp3CountM\t" +
                    "bp1CountEQ\tbp2CountEQ\tbp3CountEQ\tnCount\tdCount\t" +
                    "bp1PropM\tbp2PropM\tbp3PropM\ttotalPropM\t" +
                    "bp1PropEQ\tbp2PropEQ\tbp3PropEQ\ttotalPropEQ\tavgNPerSAM\tavgDPerSAM\tAvgNumAlignments\n"
        )

        for (transcriptName in samCountByTranscript.keys.sorted()) {
            val (xOrEqArray, eqArray, nCount, dCount, numAlignments) = aggregateMap[transcriptName]!!
            val xOrEqArrayBpSummary = getBpSummary(xOrEqArray)
            val eqArrayBpSummary = getBpSummary(eqArray)
            output.write(
                "${transcriptName}\t${xOrEqArray.size}\t${samCountByTranscript[transcriptName]}\t" +
                        "${xOrEqArrayBpSummary[0]}\t${xOrEqArrayBpSummary[1]}\t${xOrEqArrayBpSummary[2]}\t" +
                        "${eqArrayBpSummary[0]}\t${eqArrayBpSummary[1]}\t${eqArrayBpSummary[2]}\t${nCount}\t${dCount}\t" +
                        "${3 * (xOrEqArrayBpSummary[0].toDouble() / (samCountByTranscript[transcriptName]!! * xOrEqArray.size))}\t" +
                        "${3 * (xOrEqArrayBpSummary[1].toDouble() / (samCountByTranscript[transcriptName]!! * xOrEqArray.size))}\t" +
                        "${3 * (xOrEqArrayBpSummary[2].toDouble() / (samCountByTranscript[transcriptName]!! * xOrEqArray.size))}\t" +
                        "${(xOrEqArrayBpSummary[0] + xOrEqArrayBpSummary[1] + xOrEqArrayBpSummary[2]).toDouble() / (samCountByTranscript[transcriptName]!! * xOrEqArray.size)}\t" +
                        "${3 * (eqArrayBpSummary[0].toDouble() / (samCountByTranscript[transcriptName]!! * xOrEqArray.size))}\t" +
                        "${3 * (eqArrayBpSummary[1].toDouble() / (samCountByTranscript[transcriptName]!! * xOrEqArray.size))}\t" +
                        "${3 * (eqArrayBpSummary[2].toDouble() / (samCountByTranscript[transcriptName]!! * xOrEqArray.size))}\t" +
                        "${(eqArrayBpSummary[0] + eqArrayBpSummary[1] + eqArrayBpSummary[2]).toDouble() / (samCountByTranscript[transcriptName]!! * xOrEqArray.size)}\t" +
                        "${nCount.toDouble() / samCountByTranscript[transcriptName]!!}\t" +
                        "${dCount.toDouble() / samCountByTranscript[transcriptName]!!}\t" +
                        "${numAlignments.toDouble() / samCountByTranscript[transcriptName]!!}" +
                        "\n"
            )
        }
    }
}

/**
 * Function to count all the different base pair types for each exon.
 */
fun outputCountsFullTranscript(
    outputFile: String,
    samCountByTranscript: Map<String, Int>,
    samCountByTranscriptDenom: Int,
    aggregateMap: MutableMap<String, TranscriptBpAlignmentStats>
) {
    Utils.getBufferedWriter(outputFile).use { output ->

        output.write("TranscriptName\tExonNumber\tExonLength\tNumSams\t" +
                "num_bp_L\tnum_bp_C\tnum_bp_T\t" +
                "num_bp_L_Map\tnum_bp_L_EQ\t" +
                "prop_bp_L_Map\tprop_bp_L_EQ\t" +
                "num_bp_T_Map\tnum_bp_T_EQ\t" +
                "prop_bp_T_Map\tprop_bp_T_EQ\t" +
                "num_bp_C_codonPos1_Map\tnum_bp_C_codonPos1_EQ\t" +
                "num_bp_C_codonPos2_Map\tnum_bp_C_codonPos2_EQ\t" +
                "num_bp_C_codonPos3_Map\tnum_bp_C_codonPos3_EQ\t" +
                "prop_bp_C_codonPos1_Map\tprop_bp_C_codonPos1_EQ\t" +
                "prop_bp_C_codonPos2_Map\tprop_bp_C_codonPos2_EQ\t" +
                "prop_bp_C_codonPos3_Map\tprop_bp_C_codonPos3_EQ\t" +
                "num_total_Map\tnum_total_EQ\t" +
                "prop_total_Map\tprop_total_EQ\t" +
                "avgNPerSAM\tavgDPerSAM\tAvgNumAlignments\n")


        for (transcriptName in samCountByTranscript.keys.sorted()) {

            val (xOrEqArray, eqArray, nCount, dCount, numAlignments) = aggregateMap[transcriptName]!!



            val samCount = samCountByTranscript[transcriptName]!!
            val samCountDenom = if(samCountByTranscriptDenom != -1) samCountByTranscriptDenom else samCount
            //Parse out the exon numbers:
            val exonNumbers = parseExonBoundaries(transcriptName)

            var cdsCounter = 0
            var currentAlignmentBp = 0

            //Create total counters over the full transcript
            var totalNumML = 0
            var totalNumEqL = 0
            var totalNumMT = 0
            var totalNumEqT = 0

            var totalNumMCBp1 = 0
            var totalNumEqCBp1 = 0
            var totalNumMCBp2 = 0
            var totalNumEqCBp2 = 0
            var totalNumMCBp3 = 0
            var totalNumEqCBp3 = 0

            //Create some counters which are sums of all the different types of poitions we can see.
            var totalTranscriptLength = 0
            var totalNumL = 0
            var totalNumC = 0
            var totalNumCBp1 = 0
            var totalNumCBp2 = 0
            var totalNumCBp3 = 0
            var totalNumT = 0

            for((index,exon) in exonNumbers.withIndex()) {
                //Create temporary  counts for each different type of position within a given exon
                var numML = 0
                var numEqL = 0
                var numMT = 0
                var numEqT = 0

                var numMCBp1 = 0
                var numEqCBp1 = 0
                var numMCBp2 = 0
                var numEqCBp2 = 0
                var numMCBp3 = 0
                var numEqCBp3 = 0

                var totalLength = 0
                var numL = 0
                var numC = 0
                var numCBp1 = 0
                var numCBp2 = 0
                var numCBp3 = 0
                var numT = 0

                //Loop through each operator in the given exon string.
                //This should be 1,2 or 3 operators.
                for(exonBoundary in exon) {
                    val sizeOfExonOp = exonBoundary.length
                    totalLength+=sizeOfExonOp //Increment the total length of the exon by the size of this given exon part
                    when (exonBoundary.operator) {
                        ExonOperator.L -> {
                            //Increment the corresponding Leader sequence counts
                            numML += xOrEqArray.slice(currentAlignmentBp until currentAlignmentBp+sizeOfExonOp).sum()
                            numEqL += eqArray.slice(currentAlignmentBp until currentAlignmentBp+sizeOfExonOp).sum()
                            numL += sizeOfExonOp
                        }
                        ExonOperator.T -> {
                            //Increment the corresponding Terminator Sequence counts
                            numMT += xOrEqArray.slice(currentAlignmentBp until currentAlignmentBp+sizeOfExonOp).sum()
                            numEqT += eqArray.slice(currentAlignmentBp until currentAlignmentBp+sizeOfExonOp).sum()
                            numT += sizeOfExonOp
                        }
                        ExonOperator.C -> {
                            //For coding sequence we need to have a running position so we can count bp1, bp2 and bp3 correctly.
                            numC += sizeOfExonOp
                            for(idx in currentAlignmentBp until currentAlignmentBp+sizeOfExonOp) {
                                if(cdsCounter % 3 == 0) {
                                    //bp1
                                    numMCBp1 += xOrEqArray[idx]
                                    numEqCBp1 += eqArray[idx]
                                    numCBp1 ++
                                }
                                else if (cdsCounter % 3 == 1) {
                                    //bp2
                                    numMCBp2 += xOrEqArray[idx]
                                    numEqCBp2 += eqArray[idx]
                                    numCBp2++
                                }
                                else if(cdsCounter % 3 ==2) {
                                    //Bp3
                                    numMCBp3 += xOrEqArray[idx]
                                    numEqCBp3 += eqArray[idx]
                                    numCBp3++
                                }

                                cdsCounter++
                            }
                        }
                    }
                    //Shift up the current alignment Base pair by the size of the Exon part.
                    currentAlignmentBp += sizeOfExonOp
                }

                // Increase the total transcript counters correctly
                totalNumML += numML
                totalNumEqL += numEqL
                totalNumMT += numMT
                totalNumEqT += numEqT

                totalNumMCBp1 += numMCBp1
                totalNumEqCBp1 += numEqCBp1
                totalNumMCBp2 += numMCBp2
                totalNumEqCBp2 += numEqCBp2
                totalNumMCBp3 += numMCBp3
                totalNumEqCBp3 += numEqCBp3

                totalTranscriptLength +=totalLength
                totalNumL += numL
                totalNumC += numC
                totalNumCBp1 += numCBp1
                totalNumCBp2 += numCBp2
                totalNumCBp3 += numCBp3

                totalNumT += numT

                //Create total Map and EQ counts for this specific exon
                val thisExonTotalMap = numML + numMCBp1 + numMCBp2 + numMCBp3 + numMT
                val thisExonTotalEq = numEqL + numEqCBp1 + numEqCBp2 + numEqCBp3 + numEqT



                //Write out the output.  Each prop column is scaled by the number of sam files which had an alignment for that transcript
                //For the Coding sequence we need to multiply by 3 otherwise it will be 1/3 of what it really is.
                //For the average N, D and numAlignments, those values are shared across all exons as they are transcript level measurements.
                output.write("${transcriptName}\t${index}\t${totalLength}\t${samCount}\t" +
                        "${numL}\t${numC}\t${numT}\t" +
                        "${numML}\t${numEqL}\t" +
                        "${ if(numL== 0) 0.0 else numML.toDouble()/(numL * samCountDenom)}\t${ if(numL== 0) 0.0 else numEqL.toDouble()/(numL * samCountDenom)}\t" +
                        "${numMT}\t${numEqT}\t" +
                        "${ if(numT == 0) 0.0 else numMT.toDouble()/(numT * samCountDenom)}\t${ if(numT == 0) 0.0 else numEqT.toDouble()/(numT * samCountDenom)}\t" +
                        "${numMCBp1}\t${numEqCBp1}\t" +
                        "${numMCBp2}\t${numEqCBp2}\t" +
                        "${numMCBp3}\t${numEqCBp3}\t" +
                        "${ if(numCBp1 == 0) 0.0 else (numMCBp1.toDouble()) / (numCBp1 * samCountDenom)}\t${ if(numCBp1 == 0) 0.0 else (numEqCBp1.toDouble()) / (numCBp1 * samCountDenom)}\t" +
                        "${ if(numCBp2 == 0) 0.0 else (numMCBp2.toDouble()) / (numCBp2 * samCountDenom)}\t${ if(numCBp2 == 0) 0.0 else (numEqCBp2.toDouble()) / (numCBp2 * samCountDenom)}\t" +
                        "${ if(numCBp3 == 0) 0.0 else (numMCBp3.toDouble()) / (numCBp3 * samCountDenom)}\t${ if(numCBp3 == 0) 0.0 else (numEqCBp3.toDouble()) / (numCBp3 * samCountDenom)}\t" +
                        "${thisExonTotalMap}\t${thisExonTotalEq}\t" +
                        "${thisExonTotalMap.toDouble()/(totalLength * samCountDenom)}\t${thisExonTotalEq.toDouble()/(totalLength* samCountDenom)}\t" +
                        "${nCount.toDouble() / samCountDenom}\t" +
                        "${dCount.toDouble() / samCountDenom}\t" +
                        "${numAlignments.toDouble() / samCountDenom}" +
                        "\n")

            }


            //Sum up the total Mapping counts and the EQ counts and export the full transcript counts.
            val fullTranscriptTotalMap = totalNumML + totalNumMCBp1 + totalNumMCBp2 + totalNumMCBp3 + totalNumMT
            val fullTranscriptTotalEq = totalNumEqL + totalNumEqCBp1 + totalNumEqCBp2 + totalNumEqCBp3 + totalNumEqT


            output.write("${transcriptName}\tfull\t${totalTranscriptLength}\t${samCount}\t" +
                    "${totalNumL}\t${totalNumC}\t${totalNumT}\t" +
                    "${totalNumML}\t${totalNumEqL}\t" +
                    "${ if(totalNumL== 0) 0.0 else totalNumML.toDouble()/(totalNumL * samCountDenom)}\t${ if(totalNumL== 0) 0.0 else totalNumEqL.toDouble()/(totalNumL * samCountDenom)}\t" +
                    "${totalNumMT}\t${totalNumEqT}\t" +
                    "${ if(totalNumT== 0) 0.0 else totalNumMT.toDouble()/(totalNumT * samCountDenom)}\t${ if(totalNumT== 0) 0.0 else totalNumEqT.toDouble()/(totalNumT * samCountDenom)}\t" +
                    "${totalNumMCBp1}\t${totalNumEqCBp1}\t" +
                    "${totalNumMCBp2}\t${totalNumEqCBp2}\t" +
                    "${totalNumMCBp3}\t${totalNumEqCBp3}\t" +
                    "${ if(totalNumCBp1 == 0) 0.0 else (totalNumMCBp1.toDouble()) / (totalNumCBp1 * samCountDenom)}\t${ if(totalNumCBp1 == 0) 0.0 else (totalNumEqCBp1.toDouble()) / (totalNumCBp1 * samCountDenom)}\t" +
                    "${ if(totalNumCBp2 == 0) 0.0 else (totalNumMCBp2.toDouble()) / (totalNumCBp2 * samCountDenom)}\t${ if(totalNumCBp2 == 0) 0.0 else (totalNumEqCBp2.toDouble()) / (totalNumCBp2 * samCountDenom)}\t" +
                    "${ if(totalNumCBp3 == 0) 0.0 else (totalNumMCBp3.toDouble()) / (totalNumCBp3 * samCountDenom)}\t${ if(totalNumCBp3 == 0) 0.0 else (totalNumEqCBp3.toDouble()) / (totalNumCBp3 * samCountDenom)}\t" +
                    "${fullTranscriptTotalMap}\t${fullTranscriptTotalEq}\t" +
                    "${fullTranscriptTotalMap.toDouble()/(totalTranscriptLength * samCountDenom)}\t${fullTranscriptTotalEq.toDouble()/(totalTranscriptLength * samCountDenom)}\t" +
                    "${nCount.toDouble() / samCountDenom}\t" +
                    "${dCount.toDouble() / samCountDenom}\t" +
                    "${numAlignments.toDouble() / samCountDenom}" +
                    "\n")

        }
    }

}


# Full Run of 100 Assemblies
val fullTranscriptSamDir = "/workdir/hack/aimee_pipeline/arun_alignments"

val outputFullTranscriptFile = "/workdir/hack/aimee_pipeline/arun_alignments/g3ptTranscriptCountFile_90PercentCDSConservation.txt"

val outputCDSCounterFile = "/workdir/hack/aimee_pipeline/arun_alignments/g3ptTranscriptCountFile_90PercentCDSConservation_CDSConsMatrix.txt"

Sizeof.printMemoryUse()
countBpsInCDSSamDir(fullTranscriptSamDir, outputFullTranscriptFile, false, false, .9, outputCDSCounterFile)
print(countBpsInCDSSamDir(fullTranscriptSamDir, outputFullTranscriptFile, false, false, .9, outputCDSCounterFile))
Sizeof.printMemoryUse()






