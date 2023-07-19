# prediction.py
# !/usr/bin/env python3

"""
test sequences
author Long-Chen Shen
"""

import argparse
import sys
from train_test_api.test_model import eval_training
import torch
import torch.nn as nn
from utils import get_network, test_dataloader


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('model_path', type=str, help="model path")
    parser.add_argument('seq_file', type=str, help="sequence file")
    parser.add_argument('output', type=str, default=".", help="output path")

    args = parser.parse_args()

    args.weights = sys.argv[1]
    args.seq_file = sys.argv[2]
    args.output_path = sys.argv[3]
    args.with_cuda = False

    # args.weights = "./example/model_weights/maresnet-41-best.pth"
    args.seq_file = "./example/seq_file/pre_data.data"
    # args.output_path = "./example/result_file"
    # args.with_cuda = False

    args.net = "maresnet"

    cuda_condition = torch.cuda.is_available() and args.with_cuda
    args.device = torch.device("cuda:0" if cuda_condition else "cpu")
    args.b = 200
    args.warm = 0

    net = get_network(args)

    dna_test_loader = test_dataloader(seq_file=args.seq_file, num_workers=0, batch_size=args.b, shuffle=False)
    net.load_state_dict(torch.load(args.weights))
    net.eval()

    loss_function = nn.CrossEntropyLoss()
    softmax_output = nn.Softmax(dim=1)
    eval_training(net, dna_test_loader, softmax_output, args)

