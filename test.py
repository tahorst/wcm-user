import numpy as np

dt = 10
time = np.random.rand(10000) * dt
cells = np.exp(np.log(2) / dt * time)
step = 0
while step < 1000:
	for i, cell in enumerate(cells):
		cells[i] = cell * np.exp(np.log(2) / dt)
		if cells[i] > 2:
			cells[i] /= 2
			cells = np.append(cells, cells[i])

	print "%i: %f" % (step, np.mean(cells))
	step += 1