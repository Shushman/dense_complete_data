import numpy
import matplotlib.pyplot as plt
import argparse

EPSILON = 0.00001

def get_dispersion(n,d):

	disp = d*numpy.log(d) / numpy.power(n,1.0/d)
	return disp

def get_minimum_verts(d):

	minverts = numpy.power(2,d) * numpy.power(d,d/2.0) * numpy.power(numpy.log(d),d)
	return minverts

def get_work_vs_bound_lists(d,I):

	n_0 = get_minimum_verts(d)
	n_vb = n_0

	N = n_0 * numpy.power(2.0,I)

	eb_disp = get_dispersion(N,d)
	r_0 = 2.0*eb_disp + EPSILON
	r_eb = r_0

	r_max_eb = 1.0

	vb_work_list = list()
	eb_work_list = list()
	vb_opt_list = list()
	eb_opt_list = list()

	#Initial values - close to infty
	init_vb_work = n_vb * numpy.log(n_vb) + n_vb**2
	init_eb_work = N * numpy.log(N) + (N**2)*numpy.power(r_eb,d)

	fin_work = N*numpy.log(N) + N**2


	vb_opt_list.append(1/EPSILON)
	eb_opt_list.append(1/EPSILON)
	vb_work_list.append(init_vb_work)
	eb_work_list.append(init_eb_work)

	for i in range(1,I):

		#Update parameters
		n_vb = n_vb * 2.0
		vb_disp = get_dispersion(n_vb,d)
		vb_opt = 1.0 + (2*vb_disp)/(numpy.sqrt(d) - 2*vb_disp)
		vb_work = n_vb * numpy.log(n_vb) + n_vb**2
		vb_opt_list.append(vb_opt)
		vb_work_list.append(vb_work)

		if r_eb > r_max_eb:
			continue

		r_eb = r_eb * numpy.power(2.0,1.0/d)

		if r_eb > r_max_eb:
			eb_work_list.append(fin_work)
			eb_opt_list.append(1.0)
		else:
			eb_opt = 1.0 + (2*eb_disp)/(r_eb - 2*eb_disp)
			eb_work = N * numpy.log(N) + (N**2)*numpy.power(r_eb,d)
			eb_opt_list.append(eb_opt)
			eb_work_list.append(eb_work)


	#For last, just append total work
	
	vb_opt_list.append(1.0)
	vb_work_list.append(fin_work)
	

	from IPython import embed
	embed()

	return vb_work_list,vb_opt_list,eb_work_list,eb_opt_list


if __name__=='__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('--iters',type=int,required=True)
	parser.add_argument('--dims',type=int,required=True)

	args = parser.parse_args()

	vb_work,vb_opt,eb_work,eb_opt = get_work_vs_bound_lists(args.dims,args.iters)