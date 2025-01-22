import matplotlib.pyplot as plt
import numpy as np

def plot_posterior_distribution_two_params(data_param_1, data_param_2, burning = None, save = False, file_name = 'PosteriorSolution.png',**kwargs):

	param1_name = kwargs.pop('param1_name', 'Param 1')
	param2_name = kwargs.pop('param2_name', 'Param 2')
	num_bins = kwargs.pop('num_bins', 50)
	figsize = kwargs.pop('figsize', (10,6))
	
	if burning is not None:
		data_param_1 = data_param_1[burning:]
		data_param_2 = data_param_2[burning:]
	
	from scipy.stats import gaussian_kde
	kde_1 = gaussian_kde(data_param_1)
	kde_2 = gaussian_kde(data_param_2)
	x_vals_1 = np.linspace(min(data_param_1), max(data_param_1), 1000)
	x_vals_2 = np.linspace(min(data_param_2), max(data_param_2), 1000)
	kde_vals_1 = kde_1(x_vals_1)
	kde_vals_2 = kde_2(x_vals_2)
	
	fig, axs = plt.subplots(1, 2, figsize=(12, 5))
	counts1, bins1, _ = axs[0].hist(data_param_1, bins=num_bins, density=True,
	alpha=0.75, color='#0c5934',edgecolor = 'k', linewidth = 0.2)
	axs[0].plot(x_vals_1, kde_vals_1, color='k')
	axs[0].axvline(np.median(data_param_1), linestyle = '--', color = 'k')
	axs[0].set_title(f'Posterior Distribution of {param1_name}')
	axs[0].set_xlabel(param1_name)
	axs[0].set_ylabel('Density')
	axs[0].set_xlim(0, np.amax(data_param_1))
	axs[0].set_facecolor('#f7f7f2')
	
	# Melt init distribution
	counts2, bins2, _ = axs[1].hist(data_param_2, bins=num_bins, density=True,
	alpha=0.75, color='#9c1437',edgecolor = 'k', linewidth = 0.2)
	axs[1].axvline(np.median(data_param_2),linestyle = '--', color = 'k')
	axs[1].plot(x_vals_2, kde_vals_2, color='k')
	axs[1].set_title(f'Posterior Distribution of {param2_name}')
	axs[1].set_xlabel(param2_name)
	axs[1].set_ylabel('Density')
	axs[1].set_xlim(0, np.amax(data_param_2))
	axs[1].set_facecolor('#f7f7f2')
	
	if save == False:
		plt.show()
	else:
		plt.savefig(file_name,dpi = 300)
	