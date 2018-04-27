#To pickle properly, you need dill. Highly recommend.
import dill

#import these when pickling function in multiprocessing
def run_dill_encoded(payload):
    fun, args = dill.loads(payload)
    return fun(*args)
def apply_async(pool, fun, args):
    payload = dill.dumps((fun, args))
    return pool.apply_async(run_dill_encoded, (payload,))
