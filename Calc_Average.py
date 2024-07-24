# List of fractions from data
fractions = [
1.00, 1.00, 1.00, 1.00, 0.94, 0.94, 0.23, 0.23, 0.12, 0.12,
0.08, 0.08, 0.04, 0.04, 0.00, 0.00, 0.00, 0.00
]

# Calculate the average number of hydrogen bonds present in any given frame
average_saltbr = sum(fractions)
average_saltbr

print("Average number of salt bridges:", average_saltbr)