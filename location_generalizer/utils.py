import concurrent.futures
import pandas as pd
import numpy as np
from functools import partial

def parallel_func_wrapper_update_vlocation(df, func, num_cores, SEColumnMappings, homeInfo):
    df_split = np.array_split(df, num_cores)
    pool = concurrent.futures.ProcessPoolExecutor(max_workers = num_cores)
    apply_partial = partial(func, SEColumnMappings=SEColumnMappings, homeInfo=homeInfo)
    result = pd.concat(pool.map(apply_partial, df_split.copy()))
    return result

