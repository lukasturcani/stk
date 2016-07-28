def dedupe(iterable, seen=None):
    if seen is None:
        seen = set()        
    for x in iterable:
        if x not in seen:
            seen.add(x)
            yield x
            
def flatten(iterable, excluded_types={str}):
    for x in iterable:
        if hasattr(x, '__iter__') and type(x) not in excluded_types:          
            yield from flatten(x)
        else:
            yield x            
            

a = [[1],[1],[1],[[2],[2],['str']]]
