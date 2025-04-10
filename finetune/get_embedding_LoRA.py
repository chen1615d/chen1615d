import argparse
import random
import os
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
import scanpy as sc
from load import *
from datetime import datetime



####################################Settings#################################
parser = argparse.ArgumentParser(description='Drug_response_pre')
parser.add_argument('--task_name', type=str, default='deepcdr', help='task name')
parser.add_argument('--input_type', type=str, default='singlecell', choices=['singlecell', 'bulk'],
                    help='input type; default: singlecell')
parser.add_argument('--output_type', type=str, default='cell',
                    choices=['cell', 'gene', 'gene_batch', 'gene_expression'],
                    help='cell or gene embedding; default: cell')
parser.add_argument('--pool_type', type=str, default='all', choices=['all', 'max'],
                    help='pooling type of cell embedding; default: all only valid for output_type=cell')
parser.add_argument('--tgthighres', type=str, default='t4',
                    help='the targeted high resolution (start with t), fold change (start with f), or addition (start with a)')
parser.add_argument('--data_path', type=str, default='./', help='input data path')
parser.add_argument('--save_path', type=str, default='./', help='save path')
parser.add_argument('--pre_normalized', type=str, default='F', choices=['F', 'T', 'A'],
                    help='if normalized before input; default: False (F).')
parser.add_argument('--demo', action='store_true', default=False, help='if demo, only infer 10 samples')
parser.add_argument('--version', type=str, default='ce',
                    help='only valid for output_type=cell. For read depth enhancement, version=rde; For others, version=ce')
parser.add_argument('--model_path', type=str, default='None', help='pre-trained model path')
parser.add_argument('--ckpt_name', type=str, default='01B-resolution', help='checkpoint name')
parser.add_argument('--batch_size', type=int, default=4, help='Batch size for processing data; default: 16')
parser.add_argument('--device', type=str, default='cuda:1', choices=['cuda:0', 'cuda:1', 'cpu'],
                    help='Device to run the model on; default is cuda:0')

args = parser.parse_args()
ckpt_path = os.path.join(args.model_path, args.ckpt_name)
print(f"Loading model from: {ckpt_path}")
# Create timestamp in a format suitable for filenames
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

def main():
    # Set random seed
    random.seed(0)
    np.random.seed(0)
    torch.manual_seed(0)
    torch.cuda.manual_seed_all(0)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    device = torch.device(args.device if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # Load data 2025312 modify
    if args.data_path.endswith('npz'):
        gexpr_feature = scipy.sparse.load_npz(args.data_path).toarray()
        gexpr_feature = pd.DataFrame(gexpr_feature)
    elif args.data_path.endswith('h5ad'):
        gexpr_feature = sc.read_h5ad(args.data_path)
        gexpr_feature = pd.DataFrame(gexpr_feature.X)
    elif args.data_path.endswith('npy'):
        gexpr_feature = pd.DataFrame(np.load(args.data_path))
    else:
        gexpr_feature = pd.read_csv(args.data_path, index_col=0)

    if args.demo:
        gexpr_feature = gexpr_feature.iloc[:10, :]

    print(f"Input data shape: {gexpr_feature.shape}")

    # # Load data
    # if args.data_path[-3:] == 'npy':
    #     gexpr_feature = np.load(args.data_path)
    #     gexpr_feature = pd.DataFrame(gexpr_feature)
    # else:
    #     raise ValueError("Unsupported data format. Only .npy files are supported.")
    #
    # if args.demo:
    #     gexpr_feature = gexpr_feature.iloc[:10, :]
    #
    # print(f"Input data shape: {gexpr_feature.shape}")

    # Load model
    # ckpt_path = args.model_path if args.version == 'noversion' else './models/models.ckpt'
    # ckpt_path = args.model_path if args.version == 'noversion' else "/home/xt/DATA/lwc/scfoundation/fine_tunning/model/finetuned_model318.ckpt"
    ckpt_path = args.model_path if args.version == 'noversion' else "/home/xt/DATA/lwc/scfoundation/fine_tunning/model/finetuned_model_0410_1459.ckpt"
    print(f"***********: {ckpt_path}*************")
    key = None
    if args.output_type == 'cell':
        key = 'cell' if args.version == 'ce' else 'rde'
    elif args.output_type in ['gene', 'gene_batch', 'gene_expression']:
        key = 'gene'
    else:
        raise ValueError('Invalid output_type specified.')

    pretrainmodel, pretrainconfig = load_model_frommmf(ckpt_path, key)
    pretrainmodel.eval()

    # Prepare for batch processing
    BATCH_SIZE = args.batch_size
    geneexpemb = []

    for i in tqdm(range(0, gexpr_feature.shape[0], BATCH_SIZE)):
        batch_data = gexpr_feature.iloc[i:i + BATCH_SIZE].values
        # batch_data = torch.tensor(batch_data, dtype=torch.float32).cuda()
        batch_data = torch.tensor(batch_data, dtype=torch.float32).to(
            device)  # Transfer batch data to the selected device
        pretrainmodel = pretrainmodel.to(device)  # Transfer the model to the selected device

        with torch.no_grad():
            # Preprocess batch
            value_labels = batch_data > 0
            x, x_padding = gatherData(batch_data, value_labels, pretrainconfig['pad_token_id'])
            position_gene_ids, _ = gatherData(
                torch.arange(19266, device=batch_data.device).repeat(batch_data.shape[0], 1),
                value_labels, pretrainconfig['pad_token_id']
            )
            x = pretrainmodel.token_emb(torch.unsqueeze(x, 2).float(), output_weight=0)
            position_emb = pretrainmodel.pos_emb(position_gene_ids)
            x += position_emb

            # Forward pass
            geneemb = pretrainmodel.encoder(x, x_padding)

            # Pooling
            geneemb1 = geneemb[:, -1, :]
            geneemb2 = geneemb[:, -2, :]
            geneemb3, _ = torch.max(geneemb[:, :-2, :], dim=1)
            geneemb4 = torch.mean(geneemb[:, :-2, :], dim=1)

            if args.pool_type == 'all':
                geneembmerge = torch.cat([geneemb1, geneemb2, geneemb3, geneemb4], axis=1)
            elif args.pool_type == 'max':
                geneembmerge, _ = torch.max(geneemb, dim=1)
            else:
                raise ValueError('pool_type must be all or max')

            # Save batch embeddings
            geneexpemb.append(geneembmerge.detach().cpu().numpy())

    # Concatenate all embeddings and save
    geneexpemb = np.concatenate(geneexpemb, axis=0)
    print(f"Final embeddings shape: {geneexpemb.shape}")
    strname = os.path.join(args.save_path, f"{args.task_name}_{args.ckpt_name}_{args.input_type}_{args.output_type}_embedding_{args.tgthighres}_version.npy")
    strname = os.path.join(args.save_path,
                           f"{args.task_name}_{args.ckpt_name}_{args.input_type}_{args.output_type}_embedding_{args.tgthighres}_version_{timestamp}.npy")
    np.save(strname, geneexpemb)
    print(f"Embeddings saved to {strname}")

if __name__ == '__main__':
    main()
