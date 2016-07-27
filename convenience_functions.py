def dedupe(iterable, seen=None):
    if seen is None:
        seen = set()        
    for x in iterable:
        if x not in seen:
            seen.add(x)
            yield x
            
            
            

