import numpy
from logging import info
import pysam
import multiprocessing
import itertools 
#from pre_processing.shiftSize import parse_bed_for_f_r,parse_bam_for_f_r,parse_sam_for_f_r
from .pre_processing.fileParser import parse_file_by_strand, parse_file_pe

def remove_duplicate(fragments, max):
    for chr in fragments:
        fragment_chr = fragments[chr]
        fragment_chr.sort()
        idx_keep = numpy.zeros(fragment_chr.size)
        pos_pre = 0
        count = 1
        for idx,pos in enumerate(fragment_chr):
            if pos == pos_pre:
                count += 1
                if count > max:
                    idx_keep[idx] = count 
            else:
                count = 1
                pos_pre = pos
        fragments[chr] = fragment_chr[numpy.where(idx_keep==0)]
    return fragments 
        

def read_file_to_array(filename, parameter):
    ''' read file into arrays that can be used on significance testing. '''
    info ("processing %s", filename)
    move_size = parameter.window_size/2
    data_dict = {}
    for chr in parameter.chr_info:
        row_num = int(parameter.chr_info[chr]/move_size) - 1
        if row_num < 1: # if the chromosome length is less than window size. 
            row_num = 1 
        data_dict[chr] = numpy.zeros(row_num, dtype=numpy.float64)
    if parameter.file_format.endswith('pe'):
        fragments_list = [parse_file_pe[parameter.file_format](filename,parameter.input_directory)]
    else:
        forward, reverse = parse_file_by_strand[parameter.file_format](filename,parameter.input_directory)
        for chr in forward:
            forward[chr] = numpy.array(forward[chr]) + parameter.shift_dict[filename] 
        for chr in reverse:
            reverse[chr] = numpy.array(reverse[chr]) - parameter.shift_dict[filename] 
        fragments_list = [forward, reverse]
 
    for fragments in fragments_list: 
        if parameter.keep_max_dup > 0:
            fragments = remove_duplicate(fragments,parameter.keep_max_dup)
        for chr in fragments:
            for pos in fragments[chr]:
                try:
                    data_dict[chr][int(pos/move_size)] += 1
                # index out of range at the end of chr.
                except (IndexError, KeyError) as e: pass 
        
    for chr in data_dict: 
        data_dict[chr] = data_dict[chr] + numpy.roll(data_dict[chr],-1)
        data_dict[chr] = data_dict[chr] * parameter.normalization_dict[filename]
    return data_dict 
     
def read_file_to_array_wrapper(args):
    try: 
        return read_file_to_array(*args)
    except KeyboardInterrupt:
        pass
     
def read_files_to_arrays(parameter):
    
    if parameter.num_procs <2:
        for filename in parameter.get_uniq_filenames():
            parameter.array_dict[filename] = read_file_to_array(filename, parameter)
    else:
        pool = multiprocessing.Pool(parameter.num_procs)
#        p = pool.map_async(read_file_to_array_wrapper, zip(parameter.get_uniq_filenames(), itertools.repeat(parameter)),1)
        p = pool.map_async(read_file_to_array_wrapper, zip(parameter.get_uniq_filenames(), itertools.repeat(parameter)),1)
        try: results = p.get()
        except KeyboardInterrupt:
            exit(1)
        for filename, result in zip(parameter.get_uniq_filenames(), results):
            parameter.array_dict[filename] = result
            
    return 
        
def prepare_data_peak_calling(parameter):
    info ("Begin peak-calling")
    data_dict = {}
    files_remain_in_queue = parameter.get_filenames()
    for filename in parameter.chip1:
        for chr in parameter.chr_info:
            try: 
                data_dict[chr] = numpy.column_stack((data_dict[chr], parameter.array_dict[filename][chr]))
            except KeyError:
                data_dict[chr] = parameter.array_dict[filename][chr]
            if filename not in files_remain_in_queue: 
                del parameter.array_dict[filename][chr] 
        files_remain_in_queue.remove(filename)
    for filename in parameter.input1:
        for chr in parameter.chr_info:
            data_dict[chr] = numpy.column_stack((data_dict[chr], parameter.array_dict[filename][chr]))
            if filename not in files_remain_in_queue:
                del parameter.array_dict[filename][chr]        
        files_remain_in_queue.remove(filename)
    return data_dict
    
    
def prepare_data_diff_binding(parameter):
    info ("Begin differential binding analysis")
    data_dict = {}
    files_remain_in_queue = parameter.get_filenames()
    chip1_array = {}
    for filename in parameter.chip1:
        for chr in parameter.chr_info:
            try: 
                chip1_array[chr] = numpy.column_stack((chip1_array[chr], parameter.array_dict[filename][chr]))
            except KeyError:
                chip1_array[chr] = parameter.array_dict[filename][chr]
            # delete the original array. 
            if filename not in files_remain_in_queue:
                del parameter.array_dict[filename][chr]
        files_remain_in_queue.remove(filename)
    if len(parameter.input1)> 0:
        input1_array = {}
        if parameter.chip1_matched_input is True: 
            for filename in parameter.input1:
                for chr in parameter.chr_info:
                    try: 
                        input1_array[chr] = numpy.column_stack((input1_array[chr], parameter.array_dict[filename][chr]))
                    except KeyError:
                        input1_array[chr] = parameter.array_dict[filename][chr]     
                    if filename not in files_remain_in_queue:
                        del parameter.array_dict[filename][chr]
                files_remain_in_queue.remove(filename)
        else: 
            for filename in parameter.input1:
                for chr in parameter.chr_info:
                    try: 
                        input1_array[chr] += parameter.array_dict[filename][chr]
                    except KeyError:
                        input1_array[chr] = parameter.array_dict[filename][chr] 
                    if filename not in files_remain_in_queue:
                        del parameter.array_dict[filename][chr]
                files_remain_in_queue.remove(filename)

            for chr in parameter.chr_info:
                input1_array[chr] /= len(parameter.input1)
                input1_array[chr].resize(len(input1_array[chr]),1)
        for chr in parameter.chr_info:        
            chip1_array[chr] -= input1_array[chr]
            chip1_array[chr][chip1_array[chr]<0] = 0 
    
    # do the same thing for chip2 
    chip2_array = {}
    for filename in parameter.chip2:
        for chr in parameter.chr_info:
            try: 
                chip2_array[chr] = numpy.column_stack((chip2_array[chr], parameter.array_dict[filename][chr]))
            except KeyError:
                chip2_array[chr] = parameter.array_dict[filename][chr]
            if filename not in files_remain_in_queue:
                del parameter.array_dict[filename][chr]    
        files_remain_in_queue.remove(filename)
    if len(parameter.input2)> 0:
        input2_array = {}
        if parameter.chip2_matched_input is True: 
            for filename in parameter.input2:
                for chr in parameter.chr_info:
                    try: 
                        input2_array[chr] = numpy.column_stack((input2_array[chr], parameter.array_dict[filename][chr]))
                    except KeyError:
                        input2_array[chr] = parameter.array_dict[filename][chr]
                    if filename not in files_remain_in_queue:
                        del parameter.array_dict[filename][chr]
                files_remain_in_queue.remove(filename)
        else: 
            for filename in parameter.input2:
                for chr in parameter.chr_info:
                    try: 
                        input2_array[chr] += parameter.array_dict[filename][chr]
                    except KeyError:
                        input2_array[chr] = parameter.array_dict[filename][chr] 
                    if filename not in files_remain_in_queue:
                        del parameter.array_dict[filename][chr]
                files_remain_in_queue.remove(filename)
            for chr in parameter.chr_info:
                input2_array[chr] /= len(parameter.input2)
                input2_array[chr].resize(len(input2_array[chr]),1)
        for chr in parameter.chr_info:
            chip2_array[chr] -= input2_array[chr]
            chip2_array[chr][chip2_array[chr]<0] = 0    
        
    for chr in parameter.chr_info:
        data_dict[chr] =numpy.column_stack((chip1_array[chr], chip2_array[chr]))
        del chip1_array[chr]
        del chip2_array[chr]
    return data_dict
    
    
def prepare_data(parameter):
    ''' The wrapper function to prepare data. 
    Arrange the data into arrays for significance testing. '''
    if parameter.difftest == False:
        read_dict = prepare_data_peak_calling(parameter)
    else: 
        read_dict = prepare_data_diff_binding(parameter)
    # remove the array dict when it is read. 
    del parameter.array_dict 
    for chr in parameter.chr_info:
        read_array = read_dict[chr]
        read_array[numpy.where(read_array ==0)] = 1
    return read_dict
        
        
        
