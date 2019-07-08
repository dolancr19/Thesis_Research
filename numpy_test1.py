import numpy
import time

iterate=numpy.ndarray([1,4],dtype='f')
start=time.time()
for kk in range(2):
    out = numpy.array([[1,2,3,4]])
    for ii in range(10000):
        for jj in range(4):
            x=123456.789
            numpy.put(iterate,[0,jj],x)
        out=numpy.append(out,iterate, axis=0)
    numpy.save(str(kk),out)
last=time.time()
length=last-start
print(length)
