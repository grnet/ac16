def reverse(s):
    n=len(s)
    r=[-1]*n
    for i in xrange(n) : r[s[i]]=i
    return r
    
    
tr=[4,0,3,1,2]
print tr
print reverse(tr)