# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:10:44 2015

@author: yuan
"""
import math
import os
import resource
import sys
import time

class system:
	
    def __init__(self):
        pass
    
        
    def get_time(self, start_time):
        #numeric export
        #UTC time in seconds
        end_time = time.time()
        #duration by minutes
        duration = round((end_time-start_time)/60., 2)
        if duration < 0: duration = 0
        #str export
        duration_str = '{} - {} ({} min)'.format(time.ctime(start_time),time.ctime(end_time),duration)
        #time_dict
        time_dict = {'start_time':start_time, 'end_time':end_time, 
                   'duration_time':duration, 'duration':duration_str}
        return time_dict
            
    def free_space(self, pathname):
        stat = os.statvfs(pathname)
        #
        fs = round(stat.f_bfree*stat.f_bsize/(1024**3.), 2)
        print("Free space of the filesystem of {}: {} GB".format(pathname, fs))
        return fs
			
    def cpu_time(self):
        #read cpu line
        in_obj = open('/proc/stat', 'rt')
        cpu_line = in_obj.next()
        #print cpu_line
        in_obj.close()
        
        cpu_usage = {}
        #from 1-9: User time, Nice time,System time,Idle time, Waiting time
        #Hard Irq time, SoftIRQ time, Steal time
        cpu_line.rstrip("\n")
        times = cpu_line.split()
        cpu_time = math.fsum(list(map(int,times[1:])))
        cpu_usage['cpu_time'] = cpu_time
        cpu_usage['idle_time'] = int(times[4])
        #
        cpu_usage['cpu_time_perc'] = (cpu_time-cpu_usage['idle_time'])*100/cpu_time
        cpu_usage['cpu_time_perc'] = str(round(cpu_usage['cpu_time_perc'],2))+'%'
        #
        cpu_usage['user_time_perc'] = (int(times[1])+int(times[2]))*100/cpu_time
        cpu_usage['user_time_perc'] = str(round(cpu_usage['user_time_perc'],2))+'%'
        #
        cpu_usage['system_time_perc'] = (int(times[3])+int(times[6])++int(times[7]))*100/cpu_time
        cpu_usage['system_time_perc'] = str(round(cpu_usage['system_time_perc'],2))+'%'
        #print times
        #print cpu_usage
        return cpu_usage
   
    def cpu_cores(self):
        #read cpu line
        cores_num = -1
        in_obj = open('/proc/stat', 'rt')
        for line in in_obj:
            if line.find('cpu') == 0:
                cores_num += 1
                #print line
        in_obj.close()
        #print 'CPU cores:', cores_num
        return cores_num
        
    #memory usage in MB
    def used_memory(self):
        rusage_denom = 1024.
        if sys.platform == 'darwin':
            #  in OSX the output is different units ...
            rusage_denom = rusage_denom * rusage_denom
        #return by KB
        ram_used = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        ram_used = str(round(ram_used/rusage_denom,2))+' MB'
        return ram_used

    def memory_usage(self):
        #read memory line
        in_obj = open('/proc/meminfo', 'rt')
        for line in in_obj:
            line = line.rstrip("\n")
            items = line.split()
            if line.find('MemTotal') == 0:
                memTotal = int(items[1])
                #print memTotal #by kB
            elif line.find('MemFree') == 0:
                memFree = int(items[1])
                #print memFree # by kB
        in_obj.close()        
        
        #memory_usage:
        ram_usage = str(round((memTotal-memFree)*100/memTotal,2))+'%'
        memTotal = str(round(memTotal/1024.**2, 2))+' GB'
        #return percetage of used memory by % and total RAM by GB
        return ram_usage, memTotal

#calculate disk usage
    def disk_usage(self,path):
        total = os.path.getsize(path)
        #print total, path
        if os.path.isdir(path):
            for filename in os.listdir(path):
                try:
                    childpath = os.path.join(path,filename)
                    #print childpath
                    total += self.disk_usage(childpath)
                except OSError:
                    pass
        #print ('{0:<7}'.format(total), path)
        return total
        
#select parameters from STD_IN
    def select_key(self, names, label='Select data types'):
        names_exp = ['['+str(i)+']' for i in range(1, len(names)+1)]
        names_exp = '\n\t'.join([':'.join(x) for x in zip(names_exp, names)])
        while(1):
            print("{}:\n\t{}".format(label, names_exp))
            input_index = sys.stdin.readline()
            try:
                input_index = int(input_index)-1
            except ValueError:
                input_index = 0
            if 0 <= input_index < len(names):
                print(names[input_index])
                return names[input_index]

#Prints out a table using the data in `rows`, which is assumed to be a
#sequence of sequences with the 0th element being the header.
# Sample usage:
#l = [ ('title 0', 'title 1', 'title 2'), ('row 0 column 0', 'row 0 column 1 which is long', 'row 0 column 2'),
#       ('row 1 column 0', 'row 1 column 1', 'row 1 column 2')]
#print_table(l)
    def print_table(self, rows):
        # - figure out column widths
        #rows are list with nested tuples
        widths = [ len(max(list(map(str,columns)), key=len)) for columns in zip(*rows) ]
        #print widths

        # - print the separator
        print('***'.join( '*' * width for width in widths ))
        # - print the header
        header_line=[format(title, "%ds" % width) for width, title in zip(widths, rows[0])]
        print(' | '.join(header_line))
        # - print the separator
        print('-+-'.join( '-' * width for width in widths ))

        # - print the data
        for row in rows[1:]:
            line=[format(str(value), "%ds" % width) for width, value in zip(widths, row)]
            print(" | ".join(line) )
        
        # - print the separator
        print('***'.join( '*' * width for width in widths ))
#############
###end
        