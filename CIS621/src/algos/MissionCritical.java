package algos;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 * CIS 621 : Algorithms and Complexity 
 * Assignment 5
 * @author Abhishek Yenpure
 */
class Machine {

	private int cost;
	private double reliability;

	public int getCost() {
		return cost;
	}

	public void setCost(int cost) {
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

	protected static int readInput(List<Machine> machines) {
		int budget = 0;
		BufferedReader buffReader = new BufferedReader(new InputStreamReader(System.in));
		String line = null;
		try {
			int linecount = 0;
			while ((line = buffReader.readLine()) != null) {
				if (linecount == 0) {
					budget = Integer.parseInt(line);
					linecount++;
				} else if (linecount == 1) {
					linecount++;
					continue;
				} else {
					String split[] = line.split(BLANK);
					Machine machine = new Machine();
					machine.setCost(Integer.parseInt(split[0]));
					machine.setReliability(Double.parseDouble(split[1]));
					machines.add(machine);
				}
			}
			buffReader.readLine();
		} catch (IOException e) {
			System.err.println("Error occurred while reading the input, " + e.getMessage());
			System.exit(0);
		}
		return budget;
	}

	private static void initSolutionSpace(double[][] sol, int[][] count, int size, int budget) {
		for (int i = 0; i <= size; i++) {
			for (int j = 0; j <= budget; j++) {
				if (i == 0)
					sol[i][j] = 1;
				else if (j == 0)
					sol[i][j] = 0;
				else
					sol[i][j] = -1;
				count[i][j] = 0;
			}
		}
	}

	private static void printSolution(int[][] count, List<Machine> machines, int budget, double maxReliability) {
		System.out.println("Maximum Reliability : " + maxReliability);
		int machine = machines.size();
		while (machine != 0) {
			int current = count[machine][budget];
			int currentCost = machines.get(machine - 1).getCost();
			System.out.println(current + " copies of machine " + machine + " of cost " + currentCost);
			budget -= current * currentCost;
			machine--;
		}
	}

	private static double getMaxReliabilityIterative(int budget, List<Machine> machines, double[][] sol,
			int[][] count) {
		for (int i = 1; i <= machines.size(); i++) {
			Machine current = machines.get(i - 1);
			for (int j = 1; j <= budget; j++) {
				double maxReliability = 0;
				for (int k = 1; k <= j / current.getCost(); k++) {
					double currentReliability = sol[i - 1][j - k * current.getCost()]
							* (1 - Math.pow(1 - current.getReliability(), k));
					if (currentReliability > maxReliability) {
						maxReliability = currentReliability;
						sol[i][j] = maxReliability;
						count[i][j] = k;
					}
				}
			}
		}
		return sol[machines.size()][budget];
	}

	private static double getMaxReliabilityMemoized(int set, int budget, List<Machine> machines, double[][] sol,
			int[][] count) {
		if (budget < 0)
			return 0;
		else if (budget == 0 && set > 0)
			return 0;
		else if (budget >= 0 && set == 0)
			return 1;
		Machine current = machines.get(set - 1);
		double maxReliability = 0;
		for (int k = 1; k <= budget / current.getCost(); k++) {
			if (sol[set - 1][budget - k * current.getCost()] == -1)
				sol[set - 1][budget - k * current.getCost()] = getMaxReliabilityMemoized(set - 1,
						budget - k * current.getCost(), machines, sol, count);
			double currentReliability = sol[set - 1][budget - k * current.getCost()]
					* (1 - Math.pow(1 - current.getReliability(), k));
			if (currentReliability > maxReliability) {
				maxReliability = currentReliability;
				sol[set][budget] = maxReliability;
				count[set][budget] = k;
			}
		}
		return maxReliability;
	}

	private static void printMemoizationStatistics(double[][] sol, int size, int budget) {
		System.out.println("\nMemoization Statistics :");
		int total = (size + 1) * (budget + 1);
		System.out.println("Total Locations : " + total);
		int used = 0;
		for (int i = 0; i <= size; i++) {
			for (int j = 0; j <= budget; j++) {
				if (sol[i][j] != -1)
					used++;
			}
		}
		System.out.println("Locations used : " + used);
		System.out.println("Utilization % : " + (used / (double) total) * 100);
	}

	public static void main(String[] args) {
		List<Machine> machines = new ArrayList<Machine>();
		int budget = readInput(machines);
		System.out.println("Budget : " + budget + "\nNumber of Machines : " + machines.size());
		double sol[][] = new double[machines.size() + 1][budget + 1];
		int count[][] = new int[machines.size() + 1][budget + 1];

		double maxReliability;
		System.out.println("\n\nIterated version :");
		initSolutionSpace(sol, count, machines.size(), budget);
		maxReliability = getMaxReliabilityIterative(budget, machines, sol, count);
		printSolution(count, machines, budget, maxReliability);

		System.out.println("\n\nMemoized version :");
		initSolutionSpace(sol, count, machines.size(), budget);
		maxReliability = getMaxReliabilityMemoized(machines.size(), budget, machines, sol, count);
		printSolution(count, machines, budget, maxReliability);
		printMemoizationStatistics(sol, machines.size(), budget);
	}
}
