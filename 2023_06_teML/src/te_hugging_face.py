from transformers import AutoTokenizer, AutoModelForMaskedLM
import torch,sys
import pandas as pd
import numpy as np
import gc

# Import the tokenizer and the model
tokenizer = AutoTokenizer.from_pretrained("InstaDeepAI/nucleotide-transformer-2.5b-multi-species")
model = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/nucleotide-transformer-2.5b-multi-species")
device = "cuda:0" if torch.cuda.is_available() else "cpu"
model = model.to(device)

# Load in data ['gene', 'avg_fpkm', 'class', 'upSeq', 'upSeq2k']
trainSeq = pd.read_csv(filepath_or_buffer=sys.argv[1], delimiter='\t')
trainSeq.shape

# specify input sequences
inputSeq = trainSeq['teSeq']

# Run model
resEmbeddings = np.empty((0, 2560))
for i in range(0, len(inputSeq), 24):
    tokens_ids = tokenizer.batch_encode_plus(inputSeq[i:i+24], return_tensors="pt", 
                                             padding=True, truncation=True)["input_ids"]
    tokens_ids = tokens_ids.to(device)
    attention_mask = tokens_ids != tokenizer.pad_token_id
    with torch.no_grad():
        torch_outs = model(
        tokens_ids,
        attention_mask=attention_mask,
        encoder_attention_mask=attention_mask,
        output_hidden_states=True)
    embeddings = torch_outs['hidden_states'][-1].detach().cpu().numpy()
    attention_mask = tokens_ids != tokenizer.pad_token_id
    attention_mask = attention_mask.cpu().numpy()
    attention_mask = np.expand_dims(attention_mask, axis = -1)
    masked_embeddings = embeddings * attention_mask  # multiply by 0 pad tokens embeddings
    sequences_lengths = np.sum(attention_mask, axis=1)
    mean_embeddings = np.sum(masked_embeddings, axis=1) / sequences_lengths
    resEmbeddings = np.append(resEmbeddings, mean_embeddings, axis=0)
    # print progress bar
    print(f"{i+1}/{len(inputSeq)}", end="\r")

np.savetxt(fname = sys.argv[2], X=resEmbeddings, delimiter='\t')