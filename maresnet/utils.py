"""
helper function
author Long-Chen Shen
"""
import os
import sys
import re
import datetime
from torch.utils.data import DataLoader
from DNADataset import DNADataset, ToTensor, DNADataset_for_prediction, ToTensor_for_prediction
import torchvision.transforms as transforms


def get_network(args):
    """ return given network
    """
    if args.net == 'maresnet':
        from models.maresnet import multi_attention_resnet
        net = multi_attention_resnet()
    else:
        print('the network name you have entered is not supported yet')
        sys.exit()

    net = net.to(args.device)
    return net


def get_training_dataloader(path, batch_size=16, num_workers=2, shuffle=True):
    """ return training dataloader
    Args:
        path: path to dna_seq training python dataset
        batch_size: dataloader batchsize
        num_workers: dataloader num_works
        shuffle: whether to shuffle
    Returns: train_data_loader:torch dataloader object
    """
    train_DNADataset = DNADataset(path, "train", transform=transforms.Compose([ToTensor()]))
    dna_train_loader = DataLoader(train_DNADataset, shuffle=shuffle, num_workers=num_workers, batch_size=batch_size)
    return dna_train_loader


def get_valid_dataloader(path, batch_size=16, num_workers=2, shuffle=True):
    """ return valid dataloader
    Args:
        path: path to dna_seq valid python dataset
        batch_size: dataloader batchsize
        num_workers: dataloader num_works
        shuffle: whether to shuffle
    Returns: train_data_loader:torch dataloader object
    """
    valid_DNADataset = DNADataset(path, "valid", transform=transforms.Compose([ToTensor()]))
    dna_valid_loader = DataLoader(valid_DNADataset, shuffle=shuffle, num_workers=num_workers, batch_size=batch_size)
    return dna_valid_loader


def get_test_dataloader(path, batch_size=16, num_workers=2, shuffle=True):
    """ return test dataloader
    Args:
        path: path to dna_seq test python dataset
        batch_size: dataloader batchsize
        num_workers: dataloader num_works
        shuffle: whether to shuffle
    Returns: dna_training_loader:torch dataloader object
    """
    test_DNADataset = DNADataset(path, "test", transform=transforms.Compose([ToTensor()]))
    dna_test_loader = DataLoader(test_DNADataset, shuffle=shuffle, num_workers=num_workers, batch_size=batch_size)
    return dna_test_loader


def test_dataloader(seq_file, batch_size=16, num_workers=2, shuffle=False):

    test_DNADataset = DNADataset_for_prediction(seq_file, transform=transforms.Compose([ToTensor_for_prediction()]))
    dna_test_loader = DataLoader(test_DNADataset, shuffle=shuffle, num_workers=num_workers, batch_size=batch_size)
    return dna_test_loader


def most_recent_folder(net_weights, fmt):
    """
        return most recent created folder under net_weights
        if no none-empty folder were found, return empty folder
    """
    # get subfolders in net_weights
    folders = os.listdir(net_weights)

    # filter out empty folders
    folders = [f for f in folders if len(os.listdir(os.path.join(net_weights, f)))]
    if len(folders) == 0:
        return ''

    # sort folders by folder created time
    folders = sorted(folders, key=lambda f: datetime.datetime.strptime(f, fmt))
    return folders[-1]


def most_recent_weights(weights_folder):
    """
        return most recent created weights file
        if folder is empty return empty string
    """
    weight_files = os.listdir(weights_folder)
    if len(weights_folder) == 0:
        return ''

    regex_str = r'([A-Za-z0-9]+)-([0-9]+)-(regular|best)'

    # sort files by epoch
    weight_files = sorted(weight_files, key=lambda w: int(re.search(regex_str, w).groups()[1]))

    return weight_files[-1]


def last_epoch(weights_folder):
    weight_file = most_recent_weights(weights_folder)
    if not weight_file:
       raise Exception('no recent weights were found')
    resume_epoch = int(weight_file.split('-')[1])

    return resume_epoch


def best_auc_weights(weights_folder):
    """
        return the best acc .pth file in given folder, if no
        best acc weights file were found, return empty string
    """
    files = os.listdir(weights_folder)
    if len(files) == 0:
        return ''

    regex_str = r'([A-Za-z0-9]+)-([0-9]+)-(regular|best)'
    best_files = [w for w in files if re.search(regex_str, w).groups()[2] == 'best']
    if len(best_files) == 0:
        return ''

    best_files = sorted(best_files, key=lambda w: int(re.search(regex_str, w).groups()[1]))
    return best_files[-1]


# 网络参数数量
def get_parameter_number(net):
    total_num = sum(p.numel() for p in net.parameters())
    trainable_num = sum(p.numel() for p in net.parameters() if p.requires_grad)
    return {'Total': total_num, 'Trainable': trainable_num}


def save_best_result(df_path, pred_result_test):
    with open(os.path.join(df_path, "bestiter.pred"), 'w+') as f:
        for line in pred_result_test:
            f.write(str(line))
            f.write('\n')
