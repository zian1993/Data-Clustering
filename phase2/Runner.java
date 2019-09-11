package phase2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;

public class Runner {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		
		//Changing the system output to a txt file.
		//Creating a new printStream object, initialized with a new file that will contain the output of the program
		//PrintStream writer = new PrintStream(new File("out_"+args[0]));
		
		//Set the system output to the printStream object.
		//System.setOut(writer);
		
		System.out.println("K-Means Algorithm: ");
		
		//Printing out the name of the file.
		System.out.println("Dataset: " + args[0]);
		
		//Extract all the commandline arguments and pass them into
		//the kmeans object to begin the show.
		try {
			kmeans show = new kmeans(args[0], Integer.parseInt(args[1]), 
					Integer.parseInt(args[2]), Double.parseDouble(args[3]), 
					Integer.parseInt(args[4]));
		} catch (NumberFormatException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.out.println("Done!");
		
		
		//------Making a duplicate file. The following code can be uncommented out to make a duplicate
		//file with multiple copies of the 5th data point.
		
//		try {
//			coincidentCenters duplicates = new coincidentCenters(args[0]);
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
	}
}
