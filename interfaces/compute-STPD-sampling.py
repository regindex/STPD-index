#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, threading

Description = """
Tool to build the STPD-sampling of a text.
The input file cannot contain the characters 0, 1 or 2 which are
used internally by the PFP construction algorithm.
"""

dirname         =  os.path.dirname(os.path.abspath(__file__))

bigbwt_dirname  =  os.path.join(dirname, "pfp-src")
#src_dirname     =  os.path.join(dirname, "suff-set-src")

parse_exe     =  os.path.join(bigbwt_dirname, "pscan")
bwtparse_exe  =  os.path.join(bigbwt_dirname, "bwtparse")
pfpbwt_exe    =  os.path.join(bigbwt_dirname, "pfbwt")
#pfp_exe         =  os.path.join(src_dirname,    "pfp_suffixient")
#one_pass_exe    =  os.path.join(src_dirname,    "one-pass")
#linear_exe      =  os.path.join(src_dirname,    "linear-time")
#fm_exe          =  os.path.join(src_dirname,    "fm")

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input file name', type=str)
  parser.add_argument('-o', '--output', 
                      help='output files basepath (def. input)', default="", type=str)
  parser.add_argument('-w', '--wsize',  
                      help='PFP: sliding window size (def. 10)',     default=10,  type=int)
  parser.add_argument('-p', '--mod',    
                      help='PFP: hash modulus (def. 100)',           default=100, type=int)
  parser.add_argument('-t', '--threads',
                      help='PFP: number of threads (def. 1)', default=1,   type=int)
  '''
  parser.add_argument('-a', '--algorithm', 
                      help='suffixient set construction algorithm (one-pass,linear,fm,PFP)',
                      default="one-pass", type=str)
  parser.add_argument('-o', '--output', 
                      help='output files basepath (def. input)', default="", type=str)

  parser.add_argument('-w', '--wsize',  
                      help='PFP: sliding window size (def. 10)',     default=10,  type=int)
  parser.add_argument('-p', '--mod',    
                      help='PFP: hash modulus (def. 100)',           default=100, type=int)
  parser.add_argument('-t', '--threads',
                      help='PFP: number of threads (def. 1)', default=1,   type=int)
  '''
  args = parser.parse_args()

  # define output basepath
  if args.output == "":
    args.output = args.input

  logfile_name = args.input + ".stpd-sampling.log"
  # get main directory
  args.bigbwt_dir = os.path.split(sys.argv[0])[0]
  print("Sending logging messages to file:", logfile_name)

  with open(logfile_name,"a") as logfile:

    print("Computing the BWT and SA of:",args.input)
    start = time.time()

    # ---------- parsing of the input file
    command = "{exe} {file} -w {wsize} -p {modulus} -t {th}".format(
            exe = os.path.join(args.bigbwt_dir,parse_exe),
            wsize=args.wsize, modulus = args.mod, th=args.threads, file=args.input)

    command += " -s"
    print("==== Parsing. Command:", command)
    if(execute_command(command,logfile,logfile_name)!=True):
      return
    print("Elapsed time: {0:.4f}".format(time.time()-start))

    # ---------- compute the BWT of the parse
    command = "{exe} {file} -t {th}".format(
            exe = os.path.join(args.bigbwt_dir,bwtparse_exe),
            file=args.input, th=args.threads)
    command += " -s"
    print("==== Parse's BWT computation. Command:", command)
    if(execute_command(command,logfile,logfile_name)!=True):
      return
    print("Elapsed time: {0:.4f}".format(time.time()-start))

    # ---------- compute the BWT and SA of the text
    command = "{exe} {file} -w {wsize} -t {th}".format(
            exe = os.path.join(args.bigbwt_dir,pfpbwt_exe),
            wsize=args.wsize, th=args.threads, file=args.input)
    command += " -S" # computing the full suffix array
    print("==== Text's BWT and SA computation. Command:", command)
    if(execute_command(command,logfile,logfile_name)!=True):
      return
    print("Elapsed time: {0:.4f}".format(time.time()-start))

    # ---- delete temp files
    delete_PFP_files(args,logfile,logfile_name)

    elapsed_time = time.time() - start

  #print("The resulting suffixient set was sent to: "+args.output+".suff")
  print("### Elapsed time: {time} seconds".format(time=elapsed_time))

######### AUXILIARY FUNCTIONS #########

# execute command: return True is everything OK, False otherwise
def execute_command(command,logfile,logfile_name,env=None):
  try:
    subprocess.check_call(command.split(),stdout=logfile,stderr=logfile,env=env)
  except subprocess.CalledProcessError:
    print("Error executing command line:")
    print("\t"+ command)
    print("Check log file: " + logfile_name)
    return False
  return True

def execute_command_stdin(input_file,command,logfile,logfile_name,env=None):
  try:
    with open(input_file,"r") as i:
      subprocess.check_call(command.split(),stdin=i,stdout=logfile,stderr=logfile,env=env)
  except subprocess.CalledProcessError:
    print("Error executing command line:")
    print("\t"+ command)
    print("Check log file: " + logfile_name)
    return False
  return True

# delete PFP files
def delete_PFP_files(args,logfile,logfile_name):
    #if args.k==False:
    if True:

      print("==== Deleting temporary files.") # no need to show the command

      command = "rm -f {file}.parse_old {file}.parse {file}.dict " \
          "{file}.occ {file}.last {file}.bwlast {file}.ilist".format(file=args.input)
      if(execute_command(command,logfile,logfile_name)!=True):
        return
      for i in range(args.threads):
        command = "rm -f {file}.{i}.parse_old {file}.{i}.last".format(file=args.input, i=i)
        if(execute_command(command,logfile,logfile_name)!=True):
          return
      
      command = "rm -f {file}.sai {file}.bwsai".format(file=args.input);
      if(execute_command(command,logfile,logfile_name)!=True):
        return
      for i in range(args.threads):
        command = "rm -f {file}.{i}.sai".format(file=args.input, i=i)
        if(execute_command(command,logfile,logfile_name)!=True):
          return

##########################
if __name__ == '__main__':
    main()