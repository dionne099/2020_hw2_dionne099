#!/usr/bin/env python3

start = 0
end   = 99
divisor=7
nums = []

for num in range(0,100): 
    for b in range(0,num):
        if (num % divisor == 0):
            break
    else: 
        nums.append(num)
print("Printing out numbers from",start,"to",end, " not divisible by",divisor)
print(*nums, sep = "\n")
