# List of fractions from data
fractions = [
1.00, 1.00, 1.00, 1.00, 0.84, 0.84, 0.79, 0.79, 0.54, 0.54,
0.51, 0.51, 0.11, 0.11, 0.05, 0.05, 0.04, 0.04, 0.01, 0.00,
0.00, 0.00
]

# Calculate the average number of hydrogen bonds present in any given frame
average_saltbr = sum(fractions)
average_saltbr

print("Average number of salt bridges:", average_saltbr)