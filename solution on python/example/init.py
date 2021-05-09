import test 

m = [10, 25, 50, 100, 250]
file = open('test.txt', 'w') #очищаем файл
file.close()
for N in m:
    test.do(N)
