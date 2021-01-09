import os, shutil

files = os.listdir('output')
os.mkdir('output/completed')

for i in range(1, 500):
	if f'run{i}_table.xlsx' in files:
		shutil.copy(f'output/output_{i}.mat', f'output/completed/output_{i}.mat')