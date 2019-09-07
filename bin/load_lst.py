# -*- coding: utf-8 -*-
"""
Created on Mon 23 April 2018
Modified and adapted to Python 3 on Wed 17 July 2019

@author: Maxime Mougeot
@author: Jonas Karthein
@contact: maxime.mougeot@cern.ch
@contact: jonas.karthein@cern.ch
@license: MIT license
"""

import mmap
import os
import sys
import numpy as np
import pandas as pd

#-------------------------Create a dataframe containing the conversion table from the MCS6A manual--------------------------------------#

Time_Patch_Value = ['0', '5', '1', '1a', '2a', '22', '32','2','5b','Db','f3','43','c3','3']

conversion_dict = {'Data_Length' : pd.Series([2,4,4,6,6,6,6,6,8,8,8,8,8,8], index=Time_Patch_Value),
'Data_Lost_Bit' : pd.Series([np.nan ,np.nan ,np.nan ,np.nan ,np.nan ,np.nan ,47,np.nan ,63,np.nan ,47,63,np.nan ,63], index=Time_Patch_Value),
'Tag_Bits' : pd.Series([(np.nan,np.nan) ,(np.nan,np.nan) ,(np.nan,np.nan) ,(np.nan,np.nan) ,(40,47),(40,47),(np.nan,np.nan), (np.nan,np.nan), (48,62),(48,63),(48,63),(48,62),(48,63),(58,62)], index=Time_Patch_Value),
'Sweep_Counter': pd.Series([(np.nan,np.nan),(24,31),(np.nan,np.nan),(32,47),(32,39),(np.nan,np.nan),(40,46),(np.nan,np.nan),(32,47),(32,47),(40,46),(np.nan,np.nan),(np.nan,np.nan),(np.nan,np.nan)], index=Time_Patch_Value),
'Time_Bits': pd.Series([12,20,28,28,28,36,36,44,28,28,36,44,44,54], index=Time_Patch_Value),
'Max_Sweep_Length': pd.Series([0.0000004096,0.000105,0.027,0.027,0.027,6.872,6.872,1759.2,0.027,0.027,6.872,1759.2,1759.2,1801440], index=Time_Patch_Value)}

conversion_df = pd.DataFrame(conversion_dict)


def decode_binary(binary,time_patch):
    '''
    Read the binary part of the file by chunks and decode each chunk according to the format
    given in time_patch
    The length of a chunk is given by the time_patch

    Param : binary part of
    Return: nunpy array containing the converted data : tof, sweep, channel, edge, tag, fifo
    '''

    #-----------extract data from the dataframe--------------------------#
    data_length = int(conversion_df.loc[time_patch.decode('ascii'),'Data_Length'])
    nb_bits = 8*data_length #convert nb of bytes into nb of bits
    data_lost_bit = conversion_df.loc[time_patch.decode('ascii'),'Data_Lost_Bit']
    tag_bits = conversion_df.loc[time_patch.decode('ascii'),'Tag_Bits']
    sweep_counter = conversion_df.loc[time_patch.decode('ascii'),'Sweep_Counter']
    time_bits = int(conversion_df.loc[time_patch.decode('ascii'),'Time_Bits'])
    max_sweep_length = conversion_df.loc[time_patch.decode('ascii'),'Max_Sweep_Length']

    steps = len(binary[binary.tell():])/data_length
    first_it = True

    for i in range(int(steps)):
        byteword = binary.read(data_length)
        tof, sweep, channel, edge, tag, fifo = convert_bytes(byteword,nb_bits,
        data_lost_bit, tag_bits, sweep_counter, time_bits)
        if channel != 0 :#means for real data
            if first_it:
                converted_data = np.array([tof, sweep, channel, edge, tag, fifo])
                first_it = False
            else :
                converted_data = np.vstack((converted_data, np.array([tof, sweep, channel, edge, tag, fifo])))
    binary.close()
    return converted_data


def convert_bytes(bytearray,nb_bits,data_lost_bit,tag_bits,sweep_counter,time_bits):
    '''
    Perform the actual conversion of a single event from binary to integer numbers

    See pages  5-19 and 5-20 of FastCom MCS6A for more information on the bits meaning

    Param:
        bytearray : an array of nb_bits/8 bytes encapsulating the data of a single MCS6A stop event
        nb_bits : total number of bits on which the data is encoded
        data_lost_bit : Data lost bit. Indicates if the fifo was fullself. 1 bit index
        tag_bits : Tag bits. Tag info of a single stop event (see manual). Tuple containing the bit indexes.
        sweep_counter: Sweep number of a single stop event. Tuple containing the bit indexes.
        time_bits : Number of bits encoding the time of flight a single stop eventself.
        The tof seems to be given in an unbinned format with 100ps resolution (to be confirmed).

    Return:
        Decoded tof, sweep, channel, edge, tag, fifo
    '''
    bit_word = ''

    for bytes in reversed(bytearray):
        bit_word += '{0:08b}'.format(ord(chr(bytes)))

    #convert data lost bit always first index in the reversed array (last in the manual)
    if np.isnan(data_lost_bit):
        fifo = np.nan
        index_data_lost_bit = -1
    else:
        index_data_lost_bit = nb_bits-1-int(data_lost_bit)
        fifo = int(bit_word[index_data_lost_bit],2)

    #convert tag bit
    if np.isnan(tag_bits[0]):
        tag = np.nan
        index_high_tag = -1
        index_low_tag = -1
    else:
        index_high_tag = nb_bits-1-int(tag_bits[1])
        index_low_tag = nb_bits-1-int(tag_bits[0])
        tag = int(bit_word[index_high_tag:index_low_tag+1],2)

    #convert sweep number
    if np.isnan(sweep_counter[0]):
        sweep = np.nan
        index_high_sweep = -1
        index_low_sweep = -1
    else:
        index_high_sweep = nb_bits-1-int(sweep_counter[1])
        index_low_sweep = nb_bits-1-int(sweep_counter[0])
        sweep = int(bit_word[index_high_sweep:index_low_sweep+1],2)

    #convert time of flight
    index_high_tof = max(index_low_sweep,index_low_tag,index_data_lost_bit)+1
    index_low_tof = index_high_tof+time_bits
    tof = int(bit_word[index_high_tof:index_low_tof],2)

    #these are always there no matter the format
    channel = int(bit_word[index_low_tof+1:],2)
    edge = int(bit_word[index_low_tof],2)

    return tof, sweep-1, channel, edge, tag, fifo

def get_time_patch_and_binary(listfile):
    '''
    Memory map the list file and isolate the time_patch and the binary part of the data
    Param:
        listfile : input list file
    Return:
        mapped_file : memore map of the input listfile
        time_patch : string code indicating the format in which the data are written (see manual)
    '''
    mapped_file = mmap.mmap(listfile.fileno(), 0, access=mmap.ACCESS_READ)
    search_dict = {'section' : '[DATA]' , 'list_file_type' : 'time_patch'}
    #-----------------set file index to time patch code -----------------#
    pos_type_from = mapped_file.find(search_dict['list_file_type'].encode('ascii'))+len(search_dict['list_file_type'])+1
    mapped_file.seek(pos_type_from)
    time_patch = mapped_file.readline().strip('\r\n'.encode('ascii'))
    #-----------set file index to beginning of DATA-----------------------------------#
    pos_data_from = mapped_file.find(search_dict['section'].encode('ascii'))
    mapped_file.seek(pos_data_from)
    #---readline and there no matter what the file index should point to the beginning of the binary data
    mapped_file.readline()

    return mapped_file, time_patch

def process_lst(lst_files):
    '''
    Process each list file contained in lst_files

    Param:
        lst_files: input single or multiple list files
    Output: write the converted data to file as integer numbers and using the header 'tof sweep channel edge tag fifo'
    '''
    data_dict = {}
    full_info = False   # for regular application the channel, edge, tag and fifo info are constant so they don't have to be saved. In that case keep full_info = False

    for filename in lst_files:

        with open(filename,'rb') as listfile:

            binary, time_patch = get_time_patch_and_binary(listfile)

            if full_info:
                converted_data = decode_binary(binary,time_patch)    # np.array with tof, sweep, channel, edge, tag, fifo
                header_res ='tof,sweep,channel,edge,tag,fifo'
                np.savetxt('{}/{}.csv'.format(os.path.split(filename)[0],os.path.splitext(os.path.basename(filename))[0]),converted_data,
                fmt = '%i,%i,%i,%i,%f,%f', header = header_res)
            else:
                converted_data = pd.DataFrame(decode_binary(binary,time_patch)[:, [0,1]], columns=['tof', 'sweep'])     # saves only tof and sweep info
                converted_data.tof = converted_data.tof/10  # 100ps -> ns
                converted_data.sweep = converted_data.sweep.astype('int64')  # sweep is only int
                converted_data.to_csv('{}/{}.csv'.format(os.path.split(filename)[0],os.path.splitext(os.path.basename(filename))[0]), index=False)
            print('File {} loaded successfully!'.format(os.path.splitext(os.path.basename(filename))[0]))
            data_dict[os.path.splitext(os.path.basename(filename))[0]] = converted_data
    if full_info == False:
        return(pd.concat(data_dict, axis=1))    # convert dict of dataframes into one dataframe with two column name levels

if __name__ == '__main__':
    # process_lst(sys.argv[1:])
    lst_files = ['data/test_run376.lst']
    process_lst([lst_files[-1]])
    print('done')
