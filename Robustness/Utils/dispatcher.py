"""Pair of routines for parallelising function execution. 
Originally written for an earlier project, now hosted at https://gist.github.com/padraic-padraic/98396271d1158f0b63865a9c172e7710"""

import multiprocessing

__all__ = ["repeat_execution", "star_execution"]

N_PROCESSORS = multiprocessing.cpu_count()

def err(exp):
    print(exp)

def repeat_execution(n, f, args=[], kwargs={}):
    if n < N_PROCESSORS:
        pool = multiprocessing.Pool(n)
    else:
        pool = multiprocessing.Pool(N_PROCESSORS)
    manager = multiprocessing.Manager()
    results = manager.list()
    def pool_callback(res):
        results.append(res)
    [pool.apply_async(f, args, kwargs, callback=pool_callback) for i in range(n)] 
    pool.close()
    pool.join()
    return results

def star_execution(f, arg_list, kwarg_list):
    if len(arg_list) < N_PROCESSORS:
        pool = multiprocessing.Pool(len(arg_list))
    else:
        pool = multiprocessing.Pool(N_PROCESSORS)
    call_with = zip(arg_list,kwarg_list)
    results = []
    for args,kwargs in call_with:
        results.append(pool.apply_async(f, args, kwargs))
    pool.close()
    pool.join()
    return [r.get() for r in results]