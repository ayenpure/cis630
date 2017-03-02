package algos;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

class Machine {
	private double cost;
	private double reliability;

	public double getCost() {
		return cost;
	}

	public void setCost(double cost) {
		this.cost = cost;
	}

	public double getReliability() {
		return reliability;
	}

	public void setReliability(double reliability) {
		this.reliability = reliability;
	}

}

public class MissionCritical {

	private static final String BLANK = " ";

	protected static double readInputFile(List<Machine> machines) {
		double budget = 0;
		FileReader reader = null;
		try {
			reader = new FileReader("inSample.txt");
		} catch (FileNotFoundException e) {
			System.err.println("Error occurred in reading the input file, " + e.getMessage());
			System.exit(0);
		}
		BufferedReader buffReader = new BufferedReader(reader);
		String line = null;
		try {
			int linecount = 0;
			while ((line = buffReader.readLine()) != null) {
				if (linecount == 0) {
					budget = Double.parseDouble(line);
					linecount++;
				} else if (linecount == 1) {
					linecount++;
					continue;
				} else {
					String split[] = line.split(BLANK);
					Machine machine = new Machine();
					machine.setCost(Double.parseDouble(split[0]));
					machine.setReliability(Double.parseDouble(split[1]));
					machines.add(machine);
				}
			}
			buffReader.readLine();
		} catch (IOException e) {
			System.err.println("Error occurred in reading the input file, " + e.getMessage());
			System.exit(0);
		}
		return budget;
	}

	private static double getOptimumReliability(List<Machine> machines, double budget) {
		if (machines.isEmpty())
			return 1;
		Machine current = machines.get(0);
		if (budget < current.getCost())
			return 0;
		int maxInstances = (int) (budget / current.getCost());
		if(machines.size() == 1)
			return (1 - Math.pow(1 - current.getReliability(), maxInstances));
		double maxReliability = 0;
		List<Machine> remainingMachines = machines.subList(1, machines.size());
		for (int instances = 1; instances <= maxInstances; instances++) {
			if(machines.size() == 12)
				System.out.println("Adding instance " + instances + " among " + maxInstances + " instances");				
			double reliability = (1 - Math.pow(1 - current.getReliability(), instances))
					* getOptimumReliability(remainingMachines, budget - instances * current.getCost());
			if (reliability > maxReliability)
				maxReliability = reliability;
		}
		return maxReliability;
	}

	public static void main(String[] args) {
		List<Machine> machines = new ArrayList<Machine>();
		double budget = readInputFile(machines);
		double reliability = getOptimumReliability(machines, budget);
		System.out.println("Max reliability is " + reliability);
	}
}
