# transform_on_690datasets.py
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
from conf import settings
from utils import get_network, get_training_dataloader, get_valid_dataloader, get_test_dataloader, \
    most_recent_folder, best_auc_weights, get_parameter_number, save_best_result


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
    args.b = 64
    args.warm = 0
    args.lr = 0.0004
    args.resume = True  # 是否中断继续训练

    dataset_path = os.path.join("Dataset", "690_dataset")
    dataset_list = sorted(os.listdir(dataset_path))

    for dataset_ in dataset_list:

        patience = settings.PATIENCE_TRANSFORM
        net = get_network(args)
        # network parameters
        print(get_parameter_number(net))
        # data preprocessing:
        dataset = os.path.join(dataset_path, dataset_)
        dna_training_loader = get_training_dataloader(path=dataset, num_workers=0, batch_size=args.b, shuffle=True)
        dna_valid_loader = get_valid_dataloader(path=dataset, num_workers=0, batch_size=args.b, shuffle=False)
        dna_test_loader = get_test_dataloader(path=dataset, num_workers=0, batch_size=args.b, shuffle=False)

        loss_function = nn.CrossEntropyLoss()
        softmax_output = nn.Softmax(dim=1)

        optimizer = optim.SGD(params=net.parameters(), lr=args.lr, momentum=0.9, weight_decay=5e-4)
        # optimizer.param_groups[0]['initial_lr'] = 0.05
        train_scheduler = optim.lr_scheduler.MultiStepLR(optimizer, milestones=settings.MILESTONES_TRANSFORM, gamma=0.8)
        # 新的模型记录地址
        checkpoint_path_new = os.path.join(settings.CHECKPOINT_TRANSFER_PATH, args.net, settings.TIME_NOW, dataset_)
        recent_folder = ""
        if args.resume:
            recent_folder = most_recent_folder(os.path.join(settings.CHECKPOINT_PATH, args.net),
                                               fmt=settings.DATE_FORMAT)
            if not recent_folder:
                raise Exception('no recent folder were found')

            checkpoint_path = os.path.join(settings.CHECKPOINT_PATH, args.net, recent_folder)
        else:
            checkpoint_path = os.path.join(settings.CHECKPOINT_PATH, args.net, settings.TIME_NOW)

        # use tensorboard
        if not os.path.exists(settings.LOG_TRANSFER_DIR):
            os.mkdir(settings.LOG_TRANSFER_DIR)

        #  record the epoch
        df_path = os.path.join(settings.LOG_TRANSFER_DIR, args.net, settings.TIME_NOW, dataset_)
        if not os.path.exists(df_path):
            os.makedirs(df_path)
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

        if not os.path.exists(checkpoint_path_new):
            os.makedirs(checkpoint_path_new)
        checkpoint_path_new = os.path.join(checkpoint_path_new, '{net}-{epoch}-{type}.pth')

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
                # net.load_state_dict(torch.load(weights_path, map_location={'cuda:1': 'cuda:0'}))
                best_auc, cur_result, pred_result_test = eval_training(net, dna_valid_loader, dna_test_loader,
                                                                       loss_function, softmax_output,
                                                                       args, train_after=False)
                # save best result
                save_best_result(df_path, pred_result_test)
                best_result_str = cur_result
                print('best valid auc is {:0.4f}'.format(best_auc))

        for epoch in range(1, settings.EPOCH_TRANSFORM + 1):
            if args.resume:
                if epoch <= resume_epoch:
                    continue

            output_interval = settings.OUTPUT_INTERVAL_TRANSFER
            log_dic = train(net, dna_training_loader, optimizer, loss_function, epoch, args,
                            output_interval)

            if epoch > args.warm:
                train_scheduler.step()

            auc_valid, cur_result, pred_result_test = eval_training(net, dna_valid_loader, dna_test_loader,
                                                                    loss_function, softmax_output, args,
                                                                    epoch=epoch, df_file=df_file,
                                                                    log_dic=log_dic, train_after=True)

            # start to save best performance model after learning rate decay to 0.01
            if best_auc < auc_valid:
                weights_path = checkpoint_path_new.format(net=args.net, epoch=epoch, type='best')
                print('saving weights file to {}'.format(weights_path))
                torch.save(net.state_dict(), weights_path)
                best_auc = auc_valid
                patience = settings.PATIENCE_TRANSFORM
                # save best result
                save_best_result(df_path, pred_result_test)
                best_result_str = cur_result
                continue

            patience -= 1
            if patience == 0:
                print("The best:", dataset_, best_result_str)
                print("The end!")
                break
