package phase2;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class coincidentCenters {

	//This class has only one job.
	//It will read a file, and output the same file again, but with duplicate data points.

	//Function to read from the file and create the duplicate data points in another file.
	public coincidentCenters(String infile) throws FileNotFoundException, IOException
	{
		//Using the bufferedreader object to read a file,
		//and using the FileReader and FileWriter object was inspired and learned from StackOverflow.
		int count =0;
		int i=0;
		
		//Make the output file name.
		String outfile = "mod_" + infile;
		
		//Start the FileWriter
		FileWriter out = new FileWriter(outfile);

		try (BufferedReader reader = new BufferedReader(new FileReader(infile))) {
			String line;
			
			while ((line = reader.readLine()) != null) {
				// process the line.

				//Check which line number is being read.
				//If it is the 5th data point, then make duplicates and save to file.
				
				//If it is the first line, don't mess with it. It's just info.
				if (count==0)
					;//Do nothing

				//Otherwise, check if its a 5th data point.
				else if ((count%5)==0)
				{
					//If it is, write this line out to the file 9 more times.
					for (i=0; i<9; i++)
						out.write(line + "\n");
				}

				//Write line out to the output file, as one normally would.
				out.write(line + "\n");
				
				//Increment count
				count++;
			}
		}
		
		//Close the writer, my friend. Don't want to have bad code.
		out.close();
	}

}
