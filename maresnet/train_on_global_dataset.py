# train_on_global_dataset.py
# !/usr/bin/env	python3
"""
train network using pytorch
author Long-Chen Shen
"""
import os
import argparse
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from train_test_api.train_api import train
from train_test_api.test_api import eval_training
from torch.utils.tensorboard import SummaryWriter

from conf import settings
from utils import get_network, get_training_dataloader, get_valid_dataloader, get_test_dataloader, \
    most_recent_folder, most_recent_weights, last_epoch, best_auc_weights, get_parameter_number

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument('-net', type=str, required=True, help='net type')
    # parser.add_argument('-gpu', action='store_true', default=False, help='use gpu or not')
    # parser.add_argument('-b', type=int, default=128, help='batch size for dataloader')
    # parser.add_argument('-warm', type=int, default=1, help='warm up training phase')
    # parser.add_argument('-lr', type=float, default=0.1, help='initial learning rate')
    # parser.add_argument('-resume', action='store_true', default=False, help='resume training')
    args = parser.parse_args()
    args.with_cuda = True
    cuda_condition = torch.cuda.is_available() and args.with_cuda

    args.net = "maresnet"
    args.device = torch.device("cuda:0" if cuda_condition else "cpu")
    args.b = 128
    args.warm = 0
    args.lr = 0.002
    args.resume = False  # 是否中断继续训练

    net = get_network(args)
    patience = settings.PATIENCE
    # network parameters
    print(get_parameter_number(net))
    # data preprocessing:
    dataset = os.path.join("Dataset", "global_dataset")

    dna_training_loader = get_training_dataloader(path=dataset, num_workers=0, batch_size=args.b, shuffle=True)
    dna_valid_loader = get_valid_dataloader(path=dataset, num_workers=0, batch_size=args.b, shuffle=False)
    dna_test_loader = get_test_dataloader(path=dataset, num_workers=0, batch_size=args.b, shuffle=False)

    loss_function = nn.CrossEntropyLoss()
    softmax_output = nn.Softmax(dim=1)

    optimizer = optim.SGD(params=net.parameters(), lr=args.lr, momentum=0.9, weight_decay=5e-4)
    train_scheduler = optim.lr_scheduler.MultiStepLR(optimizer, milestones=settings.MILESTONES, gamma=0.8)
    recent_folder = ""
    if args.resume:
        recent_folder = most_recent_folder(os.path.join(settings.CHECKPOINT_PATH, args.net), fmt=settings.DATE_FORMAT)
        if not recent_folder:
            raise Exception('no recent folder were found')

        checkpoint_path = os.path.join(settings.CHECKPOINT_PATH, args.net, recent_folder)
    else:
        checkpoint_path = os.path.join(settings.CHECKPOINT_PATH, args.net, settings.TIME_NOW)

    # use tensorboard
    if not os.path.exists(settings.LOG_DIR):
        os.mkdir(settings.LOG_DIR)

    # since tensorboard can't overwrite old values
    # so the only way is to create a new tensorboard log
    writer = SummaryWriter(log_dir=os.path.join(
        settings.LOG_DIR, args.net, settings.TIME_NOW))

    #  record the epoch
    df_path = os.path.join(settings.LOG_DIR, args.net, settings.TIME_NOW)
    if not os.path.exists(df_path):
        os.mkdir(df_path)
    df_file = os.path.join(df_path, "df_log.pickle")
    if not os.path.isfile(df_file):
        df_ = pd.DataFrame(columns=["epoch", "lr", "train_loss", "train_acc",
                                    "valid_loss", "valid_acc", "valid_auc",
                                    "test_loss", "test_acc", "test_auc"])
        df_.to_pickle(df_file)
        print("log DataFrame created!")

    # create model_weights folder to save model
    if not os.path.exists(checkpoint_path):
        os.makedirs(checkpoint_path)
    checkpoint_path = os.path.join(checkpoint_path, '{net}-{epoch}-{type}.pth')

    best_auc = 0.0
    best_result_str = ""
    resume_epoch = 0
    if args.resume:
        best_weights = best_auc_weights(os.path.join(settings.CHECKPOINT_PATH, args.net, recent_folder))
        if best_weights:
            weights_path = os.path.join(settings.CHECKPOINT_PATH, args.net, recent_folder, best_weights)
            print('found best auc weights file:{}'.format(weights_path))
            print('load best training file to test auc...')
            net.load_state_dict(torch.load(weights_path))
            best_auc, cur_result, _ = eval_training(net, dna_valid_loader, dna_test_loader,
                                                    loss_function, softmax_output,
                                                    args, train_after=False)
            best_result_str = cur_result
            print('best valid auc is {:0.4f}'.format(best_auc))

        recent_weights_file = most_recent_weights(os.path.join(settings.CHECKPOINT_PATH, args.net, recent_folder))
        if not recent_weights_file:
            raise Exception('no recent weights file were found')
        weights_path = os.path.join(settings.CHECKPOINT_PATH, args.net, recent_folder, recent_weights_file)
        print('loading weights file {} to resume training.....'.format(weights_path))
        net.load_state_dict(torch.load(weights_path))

        resume_epoch = last_epoch(os.path.join(settings.CHECKPOINT_PATH, args.net, recent_folder))

    for epoch in range(1, settings.EPOCH + 1):
        if args.resume:
            if epoch <= resume_epoch:
                continue

        output_interval = settings.OUTPUT_INTERVAL
        log_dic = train(net, dna_training_loader, optimizer, loss_function, epoch, args,
                        output_interval, is_tensorboard=True, writer=writer)

        if epoch > args.warm:
            train_scheduler.step()

        auc_valid, cur_result, _ = eval_training(net, dna_valid_loader, dna_test_loader,
                                                 loss_function, softmax_output, args,
                                                 epoch=epoch, df_file=df_file,
                                                 log_dic=log_dic, train_after=True)
        # start to save best performance model after learning rate decay to 0.01
        if best_auc < auc_valid:
            weights_path = checkpoint_path.format(net=args.net, epoch=epoch, type='best')
            print('saving weights file to {}'.format(weights_path))
            torch.save(net.state_dict(), weights_path)
            best_auc = auc_valid
            patience = settings.PATIENCE
            best_result_str = cur_result
            continue

        if not epoch % settings.SAVE_EPOCH:
            weights_path = checkpoint_path.format(net=args.net, epoch=epoch, type='regular')
            print('saving weights file to {}'.format(weights_path))
            torch.save(net.state_dict(), weights_path)
        patience -= 1
        if patience == 0:
            print("The best:", best_result_str)
            print("The end!")
            break
    writer.close()
