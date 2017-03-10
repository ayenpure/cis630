package byteland;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;

public class Solution {

	private static final String BLANK = " ";

	public static void main(String[] args) {
		BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
		int numberOfCases = 0;
		try {
			numberOfCases = Integer.parseInt(reader.readLine().trim());
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		Map<Integer, List<Integer>> roadMap = null;
		for (int i = 0; i < numberOfCases; i++) {
			roadMap = getRoadMap(reader);
			int robots = processSolution(roadMap);
			System.out.println(robots);
		}
		/*try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}*/
	}

	private static int processSolution(Map<Integer, List<Integer>> roadMap) {
		Queue<Integer> queue = new LinkedList<Integer>();
		queue.add(0);
		int robots = 1;
		while (!queue.isEmpty()) {
			List<Integer> currentLists = roadMap.get(queue.remove());
			if (currentLists.isEmpty())
				continue;
			int forwards = 0;
			for (int incidence : currentLists) {
				queue.add(incidence);
				if (!roadMap.get(incidence).isEmpty())
					forwards++;
			}
			if (forwards != 0)
				robots = robots + forwards - 1;
		}
		return robots;
	}

	private static Map<Integer, List<Integer>> getRoadMap(BufferedReader reader) {
		Map<Integer, List<Integer>> roadMap = new HashMap<Integer, List<Integer>>();
		int numberOfCities = 0;
		try {
			numberOfCities = Integer.parseInt(reader.readLine().trim());
		} catch (NumberFormatException | IOException e) {
			e.printStackTrace();
		}
		for(int i = 0; i < numberOfCities; i++) {
			List<Integer> incidents = new ArrayList<Integer>();
			roadMap.put(i, incidents);
		}
		for (int i = 0; i < numberOfCities - 1; i++) {
			String split[] = null;
			try {
				split = reader.readLine().split(BLANK);
			} catch (IOException e) {
				e.printStackTrace();
			}
			Integer roadStart = Integer.parseInt(split[0]);
			Integer roadEnd = Integer.parseInt(split[1]);
			List<Integer> incidents = roadMap.get(roadStart);
			incidents.add(roadEnd);
			roadMap.put(roadStart, incidents);
		}
		return roadMap;
	}
}