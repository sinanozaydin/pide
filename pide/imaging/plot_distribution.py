import matplotlib.pyplot as plt
import numpy as np

def plot_posterior_distribution_two_params(data_param_1, data_param_2, save = False, file_name = 'PosteriorSolution.png',**kwargs):

	param1_name = kwargs.pop('param1_name', 'Param 1')
	param2_name = kwargs.pop('param2_name', 'Param 2')
	num_bins = kwargs.pop('num_bins', 50)
	figsize = kwargs.pop('figsize', (10,6))
		
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
	axs[0].axvline(np.mean(data_param_1), linestyle = '-',color = 'k')
	axs[0].set_title(f'Posterior Distribution of {param1_name}')
	axs[0].set_xlabel(param1_name)
	axs[0].set_ylabel('Density')
	axs[0].set_xlim(0, np.amax(data_param_1))
	axs[0].set_facecolor('#f7f7f2')
	
	# Parameter 2 distribution
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
		print(f'The image is saved as: {file_name}')
		
def plot_posterior_distribution_heatmap_two_params(data_param_1, data_param_2, save = False, file_name = 'PosteriorSolution.png', **kwargs):
	
	param_1_min = kwargs.pop('param_1_min', np.amin(data_param_1))
	param_1_max = kwargs.pop('param_1_max', np.amax(data_param_1))
	param_2_min = kwargs.pop('param_2_min', np.amin(data_param_2))
	param_2_max = kwargs.pop('param_2_max', np.amax(data_param_2))
	param1_name = kwargs.pop('param1_name', 'Param 1')
	param2_name = kwargs.pop('param2_name', 'Param 2')
	colormap = kwargs.pop('colormap','viridis')

	from matplotlib.colors import LogNorm
	
	# Define the histogram bins
	
	bins_param1 = np.linspace(param_1_min, param_1_max, 300)
	bins_param2 = np.linspace(param_2_min, param_2_max, 300)
	
	# Create a 2D histogram
	heatmap, xedges, yedges = np.histogram2d(data_param_1, data_param_2, bins=(bins_param1, bins_param2))
	# heatmap_masked = np.where(heatmap == 0, np.nan, heatmap)
	heatmap[heatmap == 0] = np.min(heatmap[heatmap > 0]) * 0.1
	# Plot the heatmap
	plt.figure(figsize=(8, 6))
	vmin_ = 0.0
	vmax_ = np.amax(heatmap) / 2.0
	# vmin_ = np.percentile(heatmap, 1)
	# vmax_ = np.percentile(heatmap, 99)
	
	im = plt.imshow(
		heatmap.T, 
		origin='lower', 
		extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], 
		cmap=colormap,
		vmin = vmin_,
		vmax = vmax_,
		aspect='auto',
	)
	cbar = plt.colorbar(im, label='Density')
	im.set_clim([vmin_, vmax_])
	plt.axvline(np.median(data_param_1),linestyle = "--",color = 'r',label = "Median")
	plt.axvline(np.mean(data_param_1),linestyle = "-",color = 'r', label = "Mean")
	plt.axhline(np.median(data_param_2),linestyle = "--",color = 'r')
	plt.axhline(np.mean(data_param_2),linestyle = "-",color = 'r')
	
	plt.xlabel(param1_name)
	plt.ylabel(param2_name)
	plt.title('Sample Distribution Heatmap')
	plt.tight_layout()
	if save == False:
		plt.show()
	else:
		plt.savefig(file_name, dpi = 300)
				
def plot_misfit_acceptance_rate(misfits, acceptance_rates, save = False, file_name = 'MisfitAcceptancePlot.png', **kwargs):

	fig = plt.figure(figsize=(10,6))
	ax0 = plt.subplot(121)
	ax1 = plt.subplot(122)
	
	ax0.plot(misfits,color = 'k')
	ax0.set_facecolor('#f7f7f2')
	ax0.set_ylabel('Misfit')
	
	ax1.plot(acceptance_rates * 1e2,color = 'k')
	ax1.set_facecolor('#f7f7f2')
	ax1.set_ylabel('Acceptance Rate (%)')
	
	if save == False:
		plt.show()
	else:
		plt.savefig(file_name,dpi = 300)
	
	