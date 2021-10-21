# -*- coding: utf-8 -*-
"""
Created on Mon 23 April 2018
Modified and adapted to Python 3 on Wed 17 July 2019
Modified by Lukas.Nies@cern.ch on 21/10/2021
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

class load_lst():
    '''
    Process each list file contained in array of lst_files
    Param:
        lst_files: array of list files
    '''
    def __init__(self, file_array):
        """
        Initialize the conversion dataframe and some other varaibles
        """
        self.lst_files = file_array
        #-------------------------Create a dataframe containing the conversion table from the MCS6A manual--------------------------------------#
        Time_Patch_Value = ['0', '5', '1', '1a', '2a', '22', '32','2','5b','Db','f3','43','c3','3']
        conversion_dict = {'Data_Length' : pd.Series([2,4,4,6,6,6,6,6,8,8,8,8,8,8], index=Time_Patch_Value),
        'Data_Lost_Bit' : pd.Series([np.nan ,np.nan ,np.nan ,np.nan ,np.nan ,np.nan ,47,np.nan ,63,np.nan ,47,63,np.nan ,63], index=Time_Patch_Value),
        'Tag_Bits' : pd.Series([(np.nan,np.nan) ,(np.nan,np.nan) ,(np.nan,np.nan) ,(np.nan,np.nan) ,(40,47),(40,47),(np.nan,np.nan), (np.nan,np.nan), (48,62),(48,63),(48,63),(48,62),(48,63),(58,62)], index=Time_Patch_Value),
        'Sweep_Counter': pd.Series([(np.nan,np.nan),(24,31),(np.nan,np.nan),(32,47),(32,39),(np.nan,np.nan),(40,46),(np.nan,np.nan),(32,47),(32,47),(40,46),(np.nan,np.nan),(np.nan,np.nan),(np.nan,np.nan)], index=Time_Patch_Value),
        'Time_Bits': pd.Series([12,20,28,28,28,36,36,44,28,28,36,44,44,54], index=Time_Patch_Value),
        'Max_Sweep_Length': pd.Series([0.0000004096,0.000105,0.027,0.027,0.027,6.872,6.872,1759.2,0.027,0.027,6.872,1759.2,1759.2,1801440], index=Time_Patch_Value)}
        self.conversion_df = pd.DataFrame(conversion_dict)
        self.data_dict = {}

    def convert_bytes(self,bytearray,nb_bits,data_lost_bit,tag_bits,sweep_counter,time_bits,verbose=0):
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

        if bit_word != "000000000000000000000000000000000000000000000000" and verbose>1:
            print(f"bit_word: {bit_word}")
            print(f"index_data_lost_bit: {fifo}")
            print(f"index_high_tag: {index_high_tag}")
            print(f"index_low_tag: {index_low_tag}")
            print(f"tag: {tag}")
            print(f"index_high_sweep: {index_high_sweep}")
            print(f"index_low_sweep: {index_low_sweep}")
            print(f"sweep: {sweep}")

        #convert time of flight
        index_high_tof = max(index_low_sweep,index_low_tag,index_data_lost_bit)+1
        index_low_tof = index_high_tof+time_bits
        tof = int(bit_word[index_high_tof:index_low_tof],2)

        #these are always there no matter the format
        channel = int(bit_word[index_low_tof+1:],2)
        edge = int(bit_word[index_low_tof],2)

        # if tof != 0:
        #     print(tof, sweep-1, channel, edge, tag, fifo)

        return tof, sweep-1, channel, edge, tag, fifo

    def decode_binary(self,binary,time_patch, verbose = 0):
        '''
        Read the binary part of the file by chunks and decode each chunk according to the format
        given in time_patch
        The length of a chunk is given by the time_patch
        Param : binary part of
        Return: nunpy array containing the converted data : tof, sweep, channel, edge, tag, fifo
        '''

        #-----------extract data from the dataframe--------------------------#
        data_length = int(self.conversion_df.loc[time_patch.decode('ascii'),'Data_Length'])
        nb_bits = 8*data_length #convert nb of bytes into nb of bits
        data_lost_bit = self.conversion_df.loc[time_patch.decode('ascii'),'Data_Lost_Bit']
        tag_bits = self.conversion_df.loc[time_patch.decode('ascii'),'Tag_Bits']
        sweep_counter = self.conversion_df.loc[time_patch.decode('ascii'),'Sweep_Counter']
        time_bits = int(self.conversion_df.loc[time_patch.decode('ascii'),'Time_Bits'])
        max_sweep_length = self.conversion_df.loc[time_patch.decode('ascii'),'Max_Sweep_Length']

        steps = len(binary[binary.tell():])/data_length
        first_it = True

        if verbose>1:
            print(f"Data length: {data_length}\nN bits: {nb_bits}\nData lost bit: {data_lost_bit}\n\
                    tag bits: {tag_bits}\nsweep_counter: {sweep_counter}\ntime_bits: {time_bits}\n\
                    max sweep length: {max_sweep_length}\nsteps: {steps}\n")

        # Introduce sweep_counter_overflow: in some cases, MCS6 seems to allocated only a small amount of bits for storing the sweep number.
        # In time_patch=32 this is for example only 7 bits -> can count to 128 and then resets to 0. With the overflow counter we just
        # count how many overflows happen and just at the necessary amount of sweeps to the overall sweep number
        sweep_counter_overflow = 0
        old_sweep = 0 # for detecting when overflow happens
        # loop through all bytewords
        for i in range(int(steps)):
            if verbose>0:
                if (i%(int(steps/10))==0):
                    print(f"Step {i} of {steps}.")
            byteword = binary.read(data_length)
            tof, sweep, channel, edge, tag, fifo = self.convert_bytes(byteword,nb_bits,
                data_lost_bit, tag_bits, sweep_counter, time_bits, verbose=verbose)
            # Check whether overflow happened (for example old_sweep = 127, new sweep is 0)
            # Only do for non-zero events:
            if tof != 0: 
                if verbose>1: print(f"old_sweep: {old_sweep}")
                if old_sweep > sweep:
                    sweep_counter_overflow += 1
                if verbose>1: print(f"sweep_counter_overflow: {sweep_counter_overflow}")
                old_sweep = sweep 
                # Add overflow to the sweep number (in case sweep has 7bit int -> 2**7=128)
                sweep += sweep_counter_overflow*(2**(sweep_counter[1]-sweep_counter[0]+1))
                if verbose>1: print(f"sweep: {sweep}")
            #
            if channel != 0 :#means for real data
                if first_it:
                    converted_data = np.array([tof, sweep, channel, edge, tag, fifo])
                    first_it = False
                else :
                    converted_data = np.vstack((converted_data, np.array([tof, sweep, channel, edge, tag, fifo])))
        binary.close()
        return converted_data

    def get_time_patch_and_binary(self, listfile, verbose=False):
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

        if verbose>1:
            print(f"pos_type_from: {pos_type_from}\npos_data_from: {pos_data_from}\ntime_patch: {time_patch}")

        return mapped_file, time_patch

    def process(self,verbose=0):
        """
        Perform the processing of the lst_files 
        """
        full_info = False   # for regular application the channel, edge, tag and fifo info are constant so they don't have to be saved. In that case keep full_info = False
        for filename in self.lst_files:
            with open(filename,'rb') as listfile:

                binary, time_patch = self.get_time_patch_and_binary(listfile, verbose=verbose)

                if full_info:
                    converted_data = self.decode_binary(binary,time_patch,verbose=verbose)    # np.array with tof, sweep, channel, edge, tag, fifo
                    header_res ='tof,sweep,channel,edge,tag,fifo'
                    np.savetxt('{}/{}.csv'.format(os.path.split(filename)[0],os.path.splitext(os.path.basename(filename))[0]),converted_data,
                    fmt = '%i,%i,%i,%i,%f,%f', header = header_res)
                else:
                    converted_data = pd.DataFrame(self.decode_binary(binary,time_patch,verbose)[:, [0,1]], columns=['tof', 'sweep'])     # saves only tof and sweep info
                    converted_data.tof = converted_data.tof/10  # 100ps -> ns
                    converted_data.sweep = converted_data.sweep.astype('int64')  # sweep is only int
                    converted_data.to_csv('{}/{}.csv'.format(os.path.split(filename)[0],os.path.splitext(os.path.basename(filename))[0]), index=False)
                print('File {} loaded successfully!'.format(os.path.splitext(os.path.basename(filename))[0]))
                self.data_dict[os.path.splitext(os.path.basename(filename))[0]] = converted_data
        if full_info == False:
            return(pd.concat(self.data_dict, axis=1))    # convert dict of dataframes into one dataframe with two column name levels