a = complex(2)
b = complex(1,3)
c = complex()
print(a,b,c)

a:set(b)
print(a,b)

a:set(2)
print(a,b)

c = a + b + b
print(a,b,c)

c:peq(b)
print(a,b,c)

c = a - b
print(a,b,c)

c = a * b
print(a,b,c)

a = b * c
print(a,b,c)
print(a.r, a.i)

d = a / c
print(a,d,c)

a.r = 5
a.i = 1
print(a,b,c)
