@file:DependsOn("com.github.samtools:htsjdk:2.19.0")

import htsjdk.samtools.*
import java.io.File

// ARGUMENTS --------------------------------------------------------
val inputSamFile = args[0]
val minAlign = args[1].toDouble()
val maxNM = args[2].toDouble()
//val inputSamFile = "/home/bm646/Downloads/AN21TSTL0214_test.sam"
//val minAlign = 0.99
//val maxNM = 0.2
fun passesAlignmentFilter(currentRecord: SAMRecord, minAlignLengthProp: Double, maxNMProp: Double) : Boolean {
    check(minAlignLengthProp in 0.0 .. 1.0) {"Error, minAlignmentLengthProportion is not between 0 and 1.0"}
    check(maxNMProp in 0.0 .. 1.0) {"Error, maxNMProp is not between 0 and 1.0"}
    val alignmentLength = currentRecord.alignmentEnd - currentRecord.alignmentStart
    val secondAlignLength = currentRecord.cigar.cigarElements.filter { it.operator.consumesReadBases() }.filter { !it.operator.isClipping }.map { it.length }.sum()
    val alignReadLength = currentRecord.cigar.cigarElements.filter { it.operator.consumesReadBases() }.filter { !it.operator.isClipping }.map { it.length }.sum()
    val readLength = currentRecord.cigar.cigarElements.filter { it.operator.consumesReadBases() }.map { it.length }.sum()
    val actualReadLength = currentRecord.cigar.filter { when(it.operator) {
        CigarOperator.M, CigarOperator.I, CigarOperator.S, CigarOperator.EQ, CigarOperator.X, CigarOperator.H -> true
        else -> false
    } }.map { it.length }.sum()
    val editDist = currentRecord.getIntegerAttribute("NM")
    val alignProp = (alignReadLength.toDouble()/actualReadLength)
    val nmProp = (editDist.toDouble()/alignReadLength)
    return ((alignReadLength.toDouble()/actualReadLength) > minAlignLengthProp && (editDist.toDouble()/alignReadLength) < maxNMProp)
}
val samIterator = SamReaderFactory.makeDefault()
    .validationStringency(ValidationStringency.SILENT)
    .open(File(inputSamFile)).iterator()
val records = mutableListOf<SAMRecord>()
while(samIterator.hasNext()) {
    records.add(samIterator.next())
}
val samHeader = SamReaderFactory.makeDefault()
    .validationStringency(ValidationStringency.SILENT)
    .open(File(inputSamFile)).fileHeader
val filteredRecords = records.filter { passesAlignmentFilter(it, minAlign, maxNM) }
val samWriter = SAMFileWriterFactory().makeSAMWriter(samHeader, false, File("${inputSamFile.removeSuffix(".sam")}_filtered.sam"))
filteredRecords.forEach {
    samWriter.addAlignment(it)
}
samWriter.close()
