import os
import subprocess as sp
# from functools import wraps
# import sys
# import time
# import fire

# def timer(func):
#     @wraps(func)
#     def wrap(*args, **kwargs):
#         begin_time = time.perf_counter()
#         result = func(*args, **kwargs)
#         start_time = time.perf_counter()
#         print('func:%r args:[%r, %r] took: %2.4f sec' % (func.__name__, args, kwargs, start_time - begin_time))
#         return result
#     return wrap
# @timer
def extract_bin_long_fastqs(out_dir, threads):
    sam_name = os.path.join(out_dir, "long_read.sam")
    plasmidfile = os.path.join(out_dir, "plasmid_long.fastq")
    cmd=f"samtools  view -@ {threads} {sam_name} | awk '{{if((($3 ~ /plas/)&& ($2 == \"0\"|| $2 == \"16\"))||($2 == \"4\")) print \"@\"$1\"\\n\"$10\"\\n+\"$1\"\\n\"$11}}' > {plasmidfile}"
    # print(cmd)
    sp.run(cmd,shell=True)

# if __name__ == '__main__':
#     # Use Fire to convert function to a command line tool
#     fire.Fire(extract_bin_long_fastqs)