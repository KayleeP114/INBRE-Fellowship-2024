# List of fractions from data
fractions = [
1.00, 1.00, 0.97, 0.97, 0.82, 0.82, 0.30, 0.30, 0.19, 0.19,
0.15, 0.15, 0.15, 0.15, 0.09, 0.09, 0.04, 0.04, 0.00, 0.00,
0.00, 0.00
]

# Calculate the average number of hydrogen bonds present in any given frame
average_saltbr = sum(fractions)
average_saltbr

print("Average number of salt bridges:", average_saltbr)