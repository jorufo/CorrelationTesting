import pandas as pd
import math

def filter_on(data, column_name, value):
	try:
		return data[data[column_name] == value]
	except KeyError:
		print('Supplied column name, {}, is not present in dataset.'.format(column_name))
		exit()

def subset(data, amount):
	try:
		if type(amount) is float and amount > 0 and amount < 1:
			return data.sample(frac = amount)
		elif type(amount) is int and amount >= 1 and amount <= len(data.index):
			return data.sample(n = amount)
		else:
			raise ValueError
	except ValueError:
		print('Supplied value for `amount`, {}, is not allowed.  Parameter `amount` must be an integer or float.'.format(amount))
		if type(amount) is int:
			print('If using an integer for `amount`, the value must be greater than 0 but less than or equal to the number of rows in `data`, to return a subset of `data` with that many rows.')
		elif type(amount) is float:
			print('If using a float value for `amount`, the value must be between 0 and 1, to return that percentage of rows of `data`.')
		exit()

def generate_correlation(data, start_column=None):
	if start_column is None:
		return data.corr()
	elif type(start_column) is int and start_column >=1 and start_column <= len(data.columns):
		subset = data.drop(data.columns[1:start_column], axis = 1)
		return subset.corr()

def get_results(data):
	results = dict()
	for rowIndex, rowName in enumerate(data.index):
		for columnIndex, columnName in enumerate(data.columns):
			if columnIndex > rowIndex:
				# Ensure alphabetical order of pair of keys
				if rowName < columnName:
					results[rowName+'-'+columnName] = data.iloc[rowIndex,columnIndex]
				else:
					results[columnName+'-'+rowName] = data.iloc[rowIndex,columnIndex]
	return results

def analyze_results(results, threshold, accumulator=None):

	if accumulator == None:
		print('Please provide an empty dictionary for `accumulator`, in which to collect the results.')
		exit()

	if(threshold == None):
		print('Please provide a value for `threshold`.')
		exit()

	for key in results:
		if key not in accumulator:
			accumulator[key] = [0, list()]
		
		if results[key] > threshold:
			accumulator[key][0] += 1

		accumulator[key][1].append(results[key])

def calculate_statistics(subsets, accumulator):

	statistics = dict()

	if subsets > 0:
		for key in accumulator:
			statistics[key] = list()
			
			# Append the mean
			mean = sum(accumulator[key][1])/subsets
			statistics[key].append(mean)

			#Append the standard deviation
			valsMinusMeanSquared = 0
			for val in accumulator[key][1]:
				valsMinusMeanSquared += (val - mean)**2
			variance = valsMinusMeanSquared / subsets
			sd = math.sqrt(variance)
			statistics[key].append(sd)

	return statistics

def export_results_to_csv(results, filename):

	output_filename_with_extension = '.'.join([filename, 'csv']) 

	print(output_filename_with_extension)
	output_file = open(output_filename_with_extension, 'w')
	sorted_results = sorted(results.items(), key=lambda x: (-x[1][0], x[0]))
	for entry in sorted_results:
		output_file.write(str(entry[0]) + ',' + str(entry[1][0]) + ',' + str(entry[1][1]) + ',' + str(entry[1][2]) + '\n')
	output_file.close()


def export_permutation_results_to_csv(results, filename, data_type):
	output_filename_with_extension = '.'.join([filename, data_type, 'csv']) 

	print(output_filename_with_extension)
	output_file = open(output_filename_with_extension, 'w')
	sorted_results = sorted(results.items())
	for entry in sorted_results:
		line = ''
		count = 0
		while count < len(results[entry[0]]):
			line += str(results[entry[0]][count]) + ','
			count += 1
		line = line.rstrip(line[-1])
		output_file.write(str(entry[0]) + ',' + line + '\n')
	output_file.close()



def correlation_analysis(
	filename,
	filter_column,
	filters_to_thresholds,
	subsets,
	sample_prcnt_or_num,
	start_column,
	output_location,
	permutations = None):
	
	'''
	correlation_analysis function.  The function used to run the correlation analysis.
	:param filename: The full path and filename of the file containing the data that needs to be read and analyzed.
	:type filename: A CSV file
	:param filter_column: The column name in the data for which we want to filter.  For each filter present, the program will run that many times.
	:type filter_column: string
	:param filters_to_thresholds: A mapping for each filter, a valid entry in the filter_column column, mapped to a threshold value for that particular group.
	:type filters_to_thresholds: dictionary
	:param subsets: The number of times to generate 'subsets' subsets for each iteration of the experiment.
	:type subsets: string
	:param sample_prcnt_or_num: A value indicating either the percentage of data we want to use from the filtered list or the number of rows of data.  If a non-integer value < 1 and > 0 is supplied, then a percentage of the origanal data will be randomly sampled. If an integer is supplied, then a particular number of rows will be randomly sampled from the dataset.
	:type sample_prcnt_or_num: integer or double
	:param start_column: The column index where the data, capable of being calculated as a correlation, begins.  Used to ignore initial columns of data present before the data of interest.
	:type start_column: integer
	:param output_location: The folder on your computer where you would like the results to be output.
	:type output_location: string, fully qualified folder path
	:param permutations: (optional) Instructs the function to run in permutation mode, the number of times based on the argument provided.
	:type integer:
	'''

	#Check the validity of the permutation variable
	if permutations is not None:
		if permutations < 1 or type(permutations) is not int:
			print('`permutations` must be greater than 1 if a value is supplied.')
			exit()

	data = pd.read_csv(filename)

	for fltr in filters_to_thresholds:
		print('Processing {} data'.format(fltr))

		# Building output location and filename
		output_filename = 'invalid'
		if permutations is not None:
			output_filename = '_'.join(['permutation_test', fltr, str(permutations), str(filters_to_thresholds[fltr]), str(subsets), str(sample_prcnt_or_num)])
		else:
			output_filename = '_'.join([fltr, str(filters_to_thresholds[fltr]), str(subsets), str(sample_prcnt_or_num)])

		output_file = output_location + output_filename


		if permutations is not None:

			permutation_counts = dict()
			permutation_means = dict()
			permutation_sds = dict()
			permutation_count = 0
			while permutation_count < permutations:
				
				print('Permutation {}'.format(permutation_count))

				#Separate into the set we will permute and the set we will ignore
				to_permute = data.loc[data[filter_column].isin(list(filters_to_thresholds.keys()))].reset_index(drop=True)
				ignore = data.loc[~data[filter_column].isin(list(filters_to_thresholds.keys()))].reset_index(drop=True)

				#Get the column index for the column of interest
				column_location = data.columns.get_loc(filter_column)

				#Permute the entries for those we want to permute
				to_permute.iloc[:,column_location] = to_permute.iloc[:,column_location].sample(frac=1).values
				permuted_data = pd.concat([ignore,to_permute])

				accumulator = dict()

				filtered_data = filter_on(permuted_data, filter_column, fltr)

				count = 0
				while count < subsets:
					print('Processing subset: {}'.format(str(count+1)))
					subset_of_filter = subset(filtered_data, sample_prcnt_or_num)
					correlation = generate_correlation(subset_of_filter, start_column)
					results = get_results(correlation)
					analyze_results(results, filters_to_thresholds[fltr], accumulator)
					count += 1

				statistics = calculate_statistics(subsets, accumulator)


				#Compile all of the results into their respective data structures
				for key in accumulator:
					if key not in permutation_counts:
						permutation_counts[key] = list()
					permutation_counts[key].append(accumulator[key][0])
					
				for key in accumulator:
					if key not in permutation_means:
						permutation_means[key] = list()
					permutation_means[key].append(statistics[key][0])
					
				for key in accumulator:
					if key not in permutation_sds:
						permutation_sds[key] = list()
					permutation_sds[key].append(statistics[key][1])

				permutation_count += 1


			export_permutation_results_to_csv(permutation_counts, output_file, 'counts')
			export_permutation_results_to_csv(permutation_means, output_file, 'mean')
			export_permutation_results_to_csv(permutation_sds, output_file, 'sd')

		else:

			accumulator = dict()

			filtered_data = filter_on(data, filter_column, fltr)

			count = 0
			while count < subsets:
				print('Processing subset: {}'.format(str(count+1)))
				subset_of_filter = subset(filtered_data, sample_prcnt_or_num)
				correlation = generate_correlation(subset_of_filter, start_column)
				results = get_results(correlation)
				analyze_results(results, filters_to_thresholds[fltr], accumulator)
				count += 1

			#Get the mean and the standard deviation
			statistics = calculate_statistics(subsets, accumulator)

			gather = dict()
			for key in accumulator:
				gather[key] = list()
				gather[key].append(accumulator[key][0])
				gather[key].append(statistics[key][0])
				gather[key].append(statistics[key][1])

			export_results_to_csv(gather, output_file)



def main():

	#Non-permutation example
	correlation_analysis(
		filename = 'example_for_GitHub.csv',
		filter_column = 'group',
		filters_to_thresholds = {
			'CHR-NC' : 0.6013,
			'CHR-C' : 0.6013
		},
		subsets = 5,
		sample_prcnt_or_num = .9,
		start_column = 7,
		output_location = 'C:/some/output/location/',
		)
		

	#Permutation example
	correlation_analysis(
		filename = 'example_for_GitHub.csv',
		filter_column = 'group',
		filters_to_thresholds = {
			'CHR-NC' : 0.6013,
			'CHR-C' : 0.6013
		},
		subsets = 10,
		sample_prcnt_or_num = .9,
		start_column = 7,
		output_location = 'C:/some/output/location/',
		permutations = 3)




if __name__ == "__main__":
	main()