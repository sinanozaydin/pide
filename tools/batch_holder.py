if len(args_input) == 1:

			print('You have to run the program with options:')
			print('Current Methods:')
			print(' ')
			print(bcolors.RED + '1. -forward')

			print(' ')
			print(' ')
			print(bcolors.NC + 'To have more information for each method. Just run the program with the option only like SEEL_BATCH -forward ')
			sys.exit()

		if args_input[1] == '-forward':

			try:
				SEEL.composition_args_path = args_input[2]
				self.write_file_save_name = args_input[3]
			except IndexError:
				print('You need to put in these for -forward selection:')
				print('SEEL_BATCH -forward input_file.csv output_file.csv')
				sys.exit()
				
			self.home()
			self.read_cond_models()
			self.read_params()
			self.check_composition()

			self.calculate_density_solid()
			self.calculate_conductivity(method = 'index')

			self.write_results()

		elif args_input[1] == '-forward_multiple':

			try:
				SEEL.composition_args_path = args_input[2]
				self.pt_file = args_input[3]
				self.write_file_save_name = args_input[4]
			except IndexError:
				print('You need to put in these for -forward_multiple selection:')
				print('SEEL_BATCH -forward_multiple input_file.csv pt_file.csv output_file.csv')
				sys.exit()

			self.home()
			self.read_cond_models()
			self.read_params()
			self.check_composition()

			self.read_pt_file()
			self.dupliSEEL_composition_for_pt_file()

			self.calculate_density_solid()

			self.calculate_conductivity(method = 'array')
			self.write_results()