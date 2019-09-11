package phase1;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Random;

public class kmeans {

	//Creating the required variables.
	private String file;
	private int clusters;
	private int iterations;
	private double threshold;
	private int runs;
	private int actualrun;
	private Random generator;
	private ArrayList<data> alldata;
	private ArrayList<data> centers;
	private ArrayList<Double> sqerrors;
	private ArrayList<ArrayList<data>> clusterdata;
	private ArrayList<Integer> randnums;
	private ArrayList<Double> eachrun;

	public kmeans(String a, int b, int c, double d, int e) throws FileNotFoundException, IOException
	{
		//Setting all the variables to the incoming values.
		file = a;
		clusters = b;
		iterations = c;
		threshold = d;
		runs = e;
		actualrun = 0;

		//Initializing other variables.
		alldata= new ArrayList<data>();
		centers = new ArrayList<data>();
		clusterdata = new ArrayList<ArrayList<data>>();
		sqerrors = new ArrayList<Double>();
		randnums = new ArrayList<Integer>();
		eachrun = new ArrayList<Double>();

		//Starting the random number generator.
		generator = new Random();

		//Lets run the show my friends.
		runTheShow();
	}

	public void runTheShow() throws FileNotFoundException, IOException
	{
		//Read the file and set up data points.
		readFile();

		//Initializing the clusters
		initializeClusters();

		//Run the algorithm R number of times.
		for (actualrun = 0; actualrun<runs; actualrun++)
		{
			//Display the run information.
			System.out.println("Run #: "+ (actualrun+1));

			//Chosing the initial centers randomly.
			randomCenters();

			//Start running the algorithm.
			runLikeTheWind();

			//Print out an empty line for aesthetics
			System.out.println("");
		}

		//Printing the best run information
		System.out.println("Best run #: " + (findLowestError(eachrun)+1)+ " : SSE value: " + (eachrun.get(findLowestError(eachrun))));
	}

	//This function will run the kmeans algorithm.
	//1. Assign each point to its closest centroid.
	//2. Recompute the centroid of each cluster.
	private void runLikeTheWind()
	{
		//Some variables we might need.
		int counter = 0;
		double lastrun = Double.POSITIVE_INFINITY;
		double thisrun=0.0;

		//Run the algorithm while the set number of iterations is not reached,
		//or, while the SSE values of the iterations do not converge.
		while (counter<iterations)
		{
			//For each data point, calculate the squared error for all the centroids.
			//Then assign the data point to the cluster with the centroid that produces
			//the lowest squared error for each particular data value.
			computeErrorsAndAssign();

			//-----Graduate level code only...
			//Check if any of the clusters are empty.

			//If so, then locate the point with the highest SSE, and make this point the center of the empty cluster.

			//Retrieve the SSE value for this run.
			thisrun = calculateSSE();

			//Print the SSE value for this run to the console.
			System.out.println("Iteration "+(counter+1)+": "+thisrun);

			//Break if the SSE values have converged, i.e. last run is equal to this run.
			if (((lastrun-thisrun)/lastrun)<threshold)
				break;
			else
				lastrun = thisrun;

			//Recompute the centroid values to be used for the next run
			recomputeCentroids();

			//Increment the counter to keep track of iterations.
			counter++;
		}

		//Saving this runs best SSE value.
		eachrun.add(thisrun);
	}

	//This function will compute the squared error for each data point and centroid,
	//and will then assign each data point to the right cluster based on the lowest
	//computed squared error.
	public void computeErrorsAndAssign()
	{
		//Some variables we might need.
		int i=0;
		int z=0;
		int x=0;
		double sqerr = 0.0;
		double val = 0.0;
		ArrayList<Double> errors = new ArrayList<Double>();

		//Clear previous cluster and error data.
		clearClusters();
		sqerrors.clear();

		for (i=0; i<alldata.size();i++)
		{
			//Clear the variables.
			errors.clear();
			sqerr=0.0;

			for (z=0; z<centers.size();z++)
			{
				//Compute the squared error for each centroid and save it in an ArrayList.
				for (x=0; x<alldata.get(i).getData().size(); x++)
				{
					val = alldata.get(i).getData().get(x) - centers.get(z).getData().get(x);
					sqerr+= val*val;
				}

				//Save the sqerror for this centroid into ArrayList.
				errors.add(sqerr);
				//Clear variable for next centroid calculation
				sqerr = 0.0;
			}
			//Find the lowest calculated sqerror,
			//Save the value, 
			//And assign the datapoint to that respective cluster.
			sqerrors.add(errors.get(findLowestError(errors)));
			clusterdata.get(findLowestError(errors)).add(alldata.get(i));
		}

	}

	public double calculateSSE()
	{
		double val = 0.0;

		for (int i=0; i<sqerrors.size(); i++)
			val += sqerrors.get(i);

		return val;
	}

	//This function will recompute the values of the centroid attributes based on the 
	//data points assigned to its respective cluster.
	public void recomputeCentroids()
	{
		//For each centroid, recompute its attributes.
		for (int i=0; i<centers.size(); i++)
		{
			//Get the mean of the respective cluster.
			//Assign this mean data point as the new cluster center.
			centers.set(i, getMean(i));
		}

		//-----Graduate level code only...
		
		//Okay lets think about this.
		//First, we need another copy of the sqerrors arraylist, from where we can choose
		//the highest points, and then remove them from there too, to get the next highest.
		ArrayList<Double> temperrors = new ArrayList<Double>();
		temperrors.addAll(sqerrors);
		
		//Check if any of the clusters are empty.
		//for (int a=0; a<centers.size();a++)
		for (int a=0; a<clusterdata.size();a++)
		{
			//If so, then locate the point with the highest SSE, and make this point the center of the empty cluster.
			if (clusterdata.get(a).isEmpty())
			{
				//Locating the point in the alldata arraylist, based on the index of its sqerror in the sqerror arraylist.
				//Assigning the corresponding point as the center instead.
				centers.set(a, alldata.get(findMaxError(temperrors)));
				
				//Then remove this highest value of the squared error from our temp error arraylist, so that
				//we can now find the next highest one if we need it.
				//Assign a value of 0 to the current highest sqerror, so that in the next iteration we find the second highest error.
				temperrors.set(findMaxError(temperrors),0.0);
			}
		}

	}

	//This function will compute the mean of all the attributes of a given cluster.
	public data getMean(int a)
	{
		//Calculating the mean of a given cluster.

		//Initializing some variables.
		data val = new data();
		double temp = 0.0;
		int b=0;
		int datasize = alldata.get(0).getData().size();

		for (int i=0; i<datasize; i++)
		{
			//For each data attribute, find the mean of that data attribute
			//across all the data points of the given cluster.
			for (b=0; b<clusterdata.get(a).size(); b++)
				temp += clusterdata.get(a).get(b).getData().get(i);

			//Calculate the mean of this data attribute.
			temp = temp/clusterdata.get(a).size();

			//Save this calculated mean attribute as an attribute of the new center.
			val.add(temp);

			//Clear temp
			temp = 0.0;
		}

		//Return this newly calculated center.
		return val;
	}

	public int findLowestError(ArrayList<Double> e)
	{
		//Assign the first value as minimum temporarily.
		double min = e.get(0);
		int index = 0;

		//Find and return the index of the lowest squared error.
		for (int i=0; i<e.size(); i++)
		{
			if (e.get(i)<min)
				min = e.get(i);
		}

		return e.indexOf(min);
	}
	
	public int findMaxError(ArrayList<Double> e)
	{
		//Assign the first value as maximum temporarily.
		double max = e.get(0);
		int index = 0;

		//Find and return the index of the highest squared error.
		for (int i=0; i<e.size(); i++)
		{
			if (e.get(i)>max)
				max = e.get(i);
		}

		return e.indexOf(max);
	}


	private void initializeClusters()
	{
		//Initialize as many ArrayLists as needed, depending on the number of clusters.
		for (int i=0; i<clusters; i++)
			clusterdata.add(new ArrayList<data>());
	}

	private void clearClusters()
	{
		//Clear cluster data
		for (int i=0; i<clusters; i++)
			clusterdata.get(i).clear();
	}

	//This function will initially choose k random centers from all the data points.
	public void randomCenters()
	{

		int randnum = 0;

		//Clear previous centers data.
		centers.clear();

		//Generating random index numbers by which to pull data from the ArrayList of all data.
		//This has to be done k number of times, i.e. one for each cluster.
		for (int i=0; i<clusters; i++)
		{
			//Generate the number, ranging from 0 to the last entry in the dataset.
			randnum = generator.nextInt(alldata.size());

			//While the random number generated is contained in the Arraylist of all the generated numbers so far,
			//keep regenerating a random number until a unique one is found.
			while (randnums.contains(randnum))
				randnum = generator.nextInt(alldata.size());

			//Assign this random datapoint to the centers dataset.
			centers.add(alldata.get(randnum));
		}

		//Now we have our random centers. Lets move on.
	}

	//Function to read from the file and create the data points.
	public void readFile() throws FileNotFoundException, IOException
	{
		int count =0;

		//Using the bufferedreader object to read a file was inspired from StackOverflow.

		try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
			String line;
			data point;
			int i;
			while ((line = reader.readLine()) != null) {
				// process the line.

				//If its the first line, do nothing with it my friend. Let it slide.
				if (count==0)
					;//Do nothing

				//Otherwise:
				else
				{
					//Split the line based on the spaces,
					//Create a data point for all the attributes in each line,
					//Make sure to convert the attributes from string to double first though.
					//Save the data point to the list of data points.

					point = new data();

					for (i=0; i<line.split(" ").length;i++)
						point.add(Double.parseDouble(line.split(" ")[i]));

					alldata.add(point);
				}

				//Increment count
				count++;
			}
		}
	}
}
