package phase1;

import java.util.ArrayList;

public class data {
	
	//This class will represent each of the data points of any dataset.
	
	//Setting up the ArrayList variable, which will contain the data attributes.
	private ArrayList<Double> data;
	
	//Creating a parameterized constructor
	public data(ArrayList<Double> a)
	{
		data = a;
	}
	
	//A non-parameterized constructor, never know what you need my friend.
	public data()
	{
		data = new ArrayList<Double>();
	}
	
	public void add(double a)
	{
		data.add(a);
	}

	//Getters and setters. Might need them. Better be safe than sorry.
	public ArrayList<Double> getData() {
		return data;
	}

	public void setData(ArrayList<Double> data) {
		this.data = data;
	}
	
	

}
