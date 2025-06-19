#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, os, glob, math, platform, re, decimal

Description = """
Pipeline for testing the suffixient index locate
"""
base_path = os.path.dirname(os.path.abspath(__file__))

## BINARY FILES DECLARATION
build_sA_exe    = base_path + "/build/index/build_store_index"
locate_STPD_exe    = base_path + "/build/sources/stpd-index-src/locate"
locate_ri_exe    = base_path + "/experiments/original-r-index/build/ri-locate"
build_ri_exe     = base_path + "/competitors/original-r-index/build/ri-build"

data_folder = base_path + "/data_experiments"
r_index_filename = ".ri"
index_filenames = [[".colex-_v0","NoOpt"],[".colex-_v1","v1"],[".colex-_v2","v2"]]
pattern_lengths = [100,1000,10000,10]
repetitions = 1

# TIMEOUR DECLARATION
timeout = 86400 # 24h

# SET SOFTWARE TO GET TIME/SPACE BASED ON OS
if platform.system() == "Darwin":
    time_space_bin = "gtime"
else:
    time_space_bin = "\\usr\time\bin"

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input_file', nargs='?', help='the input file', type=str)
    parser.add_argument('--locate_STPD',  help='test locate all occurrence queries on the STPD index (def. False)',action='store_true')
    parser.add_argument('--locate_ri',  help='test locate one occurrence queries on the r-index  (def. False)',action='store_true')
    parser.add_argument('--logs_dir_name', help='Define the directory name containing the logs (no default)', type=str, required=True)
    args = parser.parse_args()

    # Specify the log directory path
    log_dir_path = os.path.join(base_path, args.logs_dir_name)
    # Check if the directory exists
    if not os.path.exists(log_dir_path):
        os.makedirs(log_dir_path) 
        print("Logs directory created successfully!")
    # Specify the logfile and result files paths
    log_file_path = os.path.join(log_dir_path, args.input_file.split('/')[-1] + ".log")
    log_res_path  = os.path.join(log_dir_path, args.input_file.split('/')[-1] + ".res")
    log_csv_path  = os.path.join(log_dir_path, args.input_file.split('/')[-1] + ".csv")

    # Open the file containing the table storing the results
    with open(log_csv_path,"w+") as res:

        # CSV file header
        header = "index,dataset,dataset size,pattern length,no. patterns,total occurrences,time,time per pattern,time per character,time per occurrence,memory peak,first phase,second phase,perc second phase\n"
        res.write(header)
        # Compute the dataset size
        command = "wc -c {file}".format(file = args.input_file)
        dataset_size = str(subprocess.check_output(command.split()).decode()).split(' ')[0]
        print("Dataset size (in bytes) =",dataset_size)

        for pattern_length in pattern_lengths:
            args.pattern_file = args.input_file + ".pat" + str(pattern_length) + ".fasta"
            print("locating patterns in:",args.pattern_file)

            command = "wc -l {file}".format(file = args.pattern_file)
            no_patterns = str(subprocess.check_output(command.split()).decode()).split(' ')[0]
            no_patterns = int(no_patterns)
            if(no_patterns%2 != 0):
                no_patterns += 1
            no_patterns = int(no_patterns/2)
            print("No patterns =",no_patterns)

            if args.locate_STPD:
                for it in index_filenames:

                    index_path = args.input_file + it[0];
                    if not os.path.exists(index_path):
                        print("File:",index_path,"not found!")

                    command = "{exe} -i {index_path} -p {pattern_file}".format(
                               exe = locate_STPD_exe, index_path = index_path, pattern_file = args.pattern_file)
                    if it[1] == "v1" or it[1] == "v2":
                        command += " -O " + it[1]
                    print("#####",command)
                    manage_file_cache(index_path)
                    manage_file_cache(args.pattern_file)

                    peak_memory = 0.0
                    tot_time = 0.0
                    tot_occs = 0.0
                    time_per_pattern = 0.0
                    time_per_character = 0.0
                    time_per_occurence = 0.0
                    time_per_first_phase = 0.0
                    time_per_second_phase = 0.0
                    second_phase_percentage = 0.0
                    for rep in range(repetitions):
                        output_str = str(subprocess.check_output(command.split())).split('\\n')[1:]
                        output_str = "".join(output_str).split(' ')

                        peak_memory += float(output_str[21])
                        tot_time += float(output_str[30])
                        tot_occs += float(output_str[37].split("Number")[0])
                        time_per_pattern += float(output_str[52])
                        time_per_character += float(output_str[58])
                        time_per_occurence += float(output_str[64])
                        time_per_first_phase += float(output_str[76])
                        time_per_second_phase += float(output_str[84])
                        second_phase_percentage += float(output_str[94].split("%")[0])

                    csv_line = it[0].split(".")[1] + "," + args.input_file + "," + str(dataset_size) + "," + str(pattern_length) + "," + str(no_patterns) + \
                               "," + str(tot_occs/repetitions) + "," + str(round(tot_time/repetitions, 2)) + "," + str(time_per_pattern/repetitions) + \
                               "," + str(time_per_character/repetitions) + "," + str(time_per_occurence/repetitions) + "," + str(peak_memory/repetitions) + \
                               "," + str(time_per_first_phase/repetitions) + "," + str(time_per_second_phase/repetitions) + "," + str(second_phase_percentage/repetitions) + "\n"
                    res.write(csv_line) 
                    res.flush()

            if args.locate_ri:

                index_path = args.input_file + r_index_filename;
                if not os.path.exists(index_path):
                    print("File:",index_path,"not found!")

                command = "{exe} {index_path} {pattern_path}".format(
                            exe = locate_ri_exe, index_path = index_path, pattern_path = args.pattern_file)
                print("#####",command)
                manage_file_cache(index_path)
                manage_file_cache(args.pattern_file)

                peak_memory = 0.0
                tot_time = 0.0
                tot_occs = 0.0
                time_per_pattern = 0.0
                time_per_character = 0.0
                time_per_occurence = 0.0
                time_per_first_phase = 0.0
                time_per_second_phase = 0.0
                second_phase_percentage = 0.0
                for rep in range(repetitions):
                    output_str = str(subprocess.check_output(command.split())).split(' ')

                    peak_memory += float(output_str[12])
                    tot_time += float(output_str[21])
                    tot_occs += float(output_str[28].split("\\n")[0])
                    time_per_pattern += float(output_str[43])
                    time_per_character += float(output_str[49])
                    time_per_occurence += float(output_str[55])
                    time_per_first_phase += float(output_str[63])
                    time_per_second_phase += float(output_str[70])
                    second_phase_percentage = float(output_str[80].split("%")[0])

                    csv_line = "ri" + "," + args.input_file + "," + str(dataset_size) + "," + str(pattern_length) + "," + str(no_patterns) + \
                               "," + str(tot_occs/repetitions) + "," + str(round(tot_time/repetitions, 2)) + "," + str(time_per_pattern/repetitions) + \
                               "," + str(time_per_character/repetitions) + "," + str(time_per_occurence/repetitions) + "," + str(peak_memory/repetitions) + \
                               "," + str(time_per_first_phase/repetitions) + "," + str(time_per_second_phase/repetitions) + "," + str(second_phase_percentage/repetitions) + "\n"
                res.write(csv_line) 
                res.flush()

def manage_file_cache(filename):
    try:
        # Open the file in read mode
        with open(filename, "r") as fd:
            # Get the length of the file
            length = os.path.getsize(filename)
            
            # Flush any buffered data to disk
            os.fdatasync(fd.fileno())
            
            # Advise the OS about dropping and loading cache
            os.posix_fadvise(fd.fileno(), 0, length, os.POSIX_FADV_DONTNEED)
            os.posix_fadvise(fd.fileno(), 0, length, os.POSIX_FADV_WILLNEED)

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    main()