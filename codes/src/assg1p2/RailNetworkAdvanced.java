package assg1p2;

import java.io.*;
import java.util.*;

public class RailNetworkAdvanced {

	private TreeMap<String,Station> stationList;
	private TreeMap<String, TrainLine> lineList;

	private TreeMap<Station, ArrayList<String>> pathFromOrigin;
	private TreeMap<Station, Integer> distanceFromOrigin;
	private TreeMap<String,TreeMap<String,Integer>> getMinDistFromOrigin = new TreeMap<>();
	   

	public RailNetworkAdvanced(String trainData, String connectionData, String linesData)
	{
		stationList = new TreeMap<>();
		lineList = new TreeMap<>();

		pathFromOrigin = new TreeMap<>();
		distanceFromOrigin = new TreeMap<>();

		try {
			readLinesData(linesData);
			readStationData(trainData);
			readConnectionData(connectionData);
		}
		catch (IOException e) {
			System.out.println("Exception encountered: " + e);
		}
	}
	
	/**
	 * Reads the CSV file containing information about the lines
	 * 
	 * @param infile
	 * @throws IOException
	 */
	public void readLinesData(String infile) throws IOException
	{
		BufferedReader in = new BufferedReader(new FileReader(infile));

		String s;
		s = in.readLine(); //skip the first line, csv headers

		while ((s = in.readLine()) != null)
		{
			String[] data = s.split(",");
			lineList.put(data[0], new TrainLine(data[0],
												data[1],
												data[2],
												data[3],
												Integer.valueOf(data[4])));
		}
		in.close();
	}
	/**
	 * Reads the CSV file containing information about the stations and 
	 * populate the TreeMap<String,Station> stationList. Each row of 
	 * the CSV file contains the name, latitude and longitude coordinates
	 * of the station.
	 * 
	 * You need to make a Station object for each row and add it to the 
	 * TreeMap<String,Station> stationList where the key value is the 
	 * name of the station (as a String).
	 * 
	 * @param infile	   the path to the file
	 * @throws IOException if the file is not found
	 */
	public void readStationData(String infile) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(infile));

		String s;
		s = in.readLine();
		String[] headers = s.split(",");

		while ((s = in.readLine()) != null)
		{
			HashMap<String, Integer> stationNumber = new HashMap<>();
			ArrayList<String> trainLine = new ArrayList<>();

			String[] data = s.split(",");

			//iterate throw train line values, if value not empty add to the list
			for (int i = 3; i < data.length - 1; i++)
			{
				if (data[i].compareTo("") != 0)
				{
					stationNumber.put(headers[i], Integer.valueOf(data[i]));
					trainLine.add(headers[i]);
					lineList.get(headers[i]).addStation(Integer.valueOf(data[i]), data[0]);
				}

			}

			stationList.put(data[0], new Station(data[0],
												Double.valueOf(data[1]),
												Double.valueOf(data[2]),
												trainLine,
												stationNumber));
		}
		in.close();
	}

	/**
	 * Reads the CSV file containing information about connectivity between 
	 * adjacent stations, and update the stations in stationList so that each
	 * Station object has a list of adjacent stations.
	 * 
	 * Each row contains two Strings separated by comma. To obtain the distance
	 * between the two stations, you need to use the latitude and longitude 
	 * coordinates together with the computeDistance() methods defined below
	 * 
	 * @param infile	   the path to the file
	 * @throws IOException if the file is not found
	 */
	public void readConnectionData(String infile) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(infile));
		String s;
		while ((s = in.readLine()) != null)
		{
			String[] data = s.split(",");
			if (!data[0].isEmpty() && !data[1].isEmpty())
			{
				//assume each pair is unique
				int distance = computeDistance(data[0], data[1]);
				stationList.get(data[0]).addNeighbour(stationList.get(data[1]), distance);
				stationList.get(data[1]).addNeighbour(stationList.get(data[0]), distance);
			}
		}
		in.close();
	}
	
	/**
	 * Given the latitude and longitude coordinates of two locations x and y, 
	 * return the distance between x and y in metres using Haversine formula,
	 * rounded down to the nearest integer.
	 * 
	 * Note that two more methods are provided below for your convenience 
	 * and you should not directly call this method
	 * 
	 * source://www.geeksforgeeks.org/haversine-formula-to-find-distance-between-two-points-on-a-sphere/
	 * 
	 * @param lat1 latitude coordinate of x
	 * @param lon1 longitude coordinate of x
	 * @param lat2 latitude coordinate of y
	 * @param lon2 longitude coordinate of y
	 * @return distance betwee
	 */
	public static int computeDistance(double lat1, double lon1, double lat2, double lon2)
	{
        // distance between latitudes and longitudes 
        double dLat = Math.toRadians(lat2 - lat1); 
        double dLon = Math.toRadians(lon2 - lon1); 
  
        // convert to radians 
        lat1 = Math.toRadians(lat1); 
        lat2 = Math.toRadians(lat2); 
  
        // apply formulae 
        double a = Math.pow(Math.sin(dLat / 2), 2) +  
                   Math.pow(Math.sin(dLon / 2), 2) *  
                   Math.cos(lat1) *  
                   Math.cos(lat2); 
        double rad = 6371.0; 
        Double c = 2 * Math.asin(Math.sqrt(a)); 
        Double distance = rad * c * 1000;
        return distance.intValue(); 
	}	
	
	/**
	 * Compute the distance between two stations in metres, where the stations
	 * are given as String objects
	 * 
	 * @param a		the first station
	 * @param b		the second station
	 * @return		the distance between the two stations in metres
	 */
	public int computeDistance(String a, String b)
	{
		Station u = stationList.get(a);
		Station v = stationList.get(b);
		return computeDistance(u.getLatitude(),u.getLongitude(),
							   v.getLatitude(),v.getLongitude());
	}
	
	/**
	 * Compute the distance between two stations in metres, where the stations
	 * are given as Station objects
	 * 
	 * @param a		the first station
	 * @param b		the second station
	 * @return		the distance between the two stations in metres
	 */
	public int computeDistance(Station a, Station b)
	{
		return computeDistance(a.getLatitude(),a.getLongitude(),
							   b.getLatitude(),b.getLongitude());
	}

	/**
	 * The method finds the shortest route (in terms of number of stops)
	 * between the origin station and the destination station.
	 * The route is returned as an ArrayList<String> containing the names of 
	 * the stations along the route, including the origin and the destination 
	 * stations.
	 * 
	 * If the route cannot be completed (there is no path between origin and
	 * destination), then return an empty ArrayList<String>
	 * 
	 * If the origin or the destination stations are not in the list of stations,
	 * return an empty ArrayList<String>. 
	 * 
	 * If the origin and the destination stations are the same, return an 
	 * ArrayList<String> containing the station.
	 * 
	 * @param origin		the starting station
	 * @param destination	the destination station
	 * @return
	 */
	public ArrayList<String> routeMinStop(String origin, String destination)
	{
		if (!stationList.containsKey(origin) || !stationList.containsKey(destination))
			return new ArrayList<>();

		//clear data sets
		pathFromOrigin.clear();

		Station originStation = stationList.get(origin);
		Station destinationStation = stationList.get(destination);

		//add original station
		pathFromOrigin.put(originStation, new ArrayList<>(List.of(origin)));

		//launch recursion, true if track been found
		if (routeMinStopRec(originStation, destinationStation, new TreeSet<>()))
			return pathFromOrigin.get(destinationStation);

		//if recursion returned false, nothing was found, return empty list
		return new ArrayList<>();
	}

	private Boolean routeMinStopRec(Station origin, Station destination, TreeSet<String> failures)
	{
		//if we're on dest already, return true
		if (origin.equals(destination))
			return true;

		//otherwise, set current origin to marked
		origin.setMarked();

		TreeMap<Station, Integer> neighbours = origin.getAdjacentStations();

		//modify distance to each neighbour
		for(Map.Entry<Station,Integer> neighbour : neighbours.entrySet())
		{
			if (!pathFromOrigin.containsKey(neighbour.getKey()) ||
					pathFromOrigin.get(neighbour.getKey()).size() > pathFromOrigin.get(origin).size() + 1)
			{
				neighbour.getKey().setUnmarked();

				ArrayList<String> list = new ArrayList<>(pathFromOrigin.get(origin));
				list.add(neighbour.getKey().getName());

				pathFromOrigin.put(neighbour.getKey(), list);
			}
		}

		Boolean result = false;

		//call each neighbour
		for(Map.Entry<Station,Integer> neighbour : neighbours.entrySet())
		{
			if (!neighbour.getKey().isMarked() && !failures.contains(neighbour.getKey().getName()))
				result |= routeMinStopRec(neighbour.getKey(), destination, failures);
		}
		return result;
	}


	/********************************STAGE 3*****************************************/
	public Station nearestStation(HashMap<Station, Integer> foundStations) {
        int smallest = Integer.MAX_VALUE;
        Station nearestStation = null; 
        for(Station s : foundStations.keySet()) {
                if(!s.isMarked() && foundStations.get(s)<smallest) {
                    smallest = foundStations.get(s);
                    nearestStation = s;
                }
                   
            }
        return nearestStation;
    }
   
   
   
    public void unMarkAll() {
        for(Station s: stationList.values())
            s.setUnmarked();
    }
   
    public void dijkstra(String origin){
        if (!stationList.containsKey(origin)) {
    		return;
        }
        
        //Initialisation
        Queue<Station> q = new LinkedList<Station>();
        TreeMap<Station, Integer> distance = new TreeMap<>();
        HashMap<Station,Integer> foundStations = new HashMap<>();
        
        //Population
        for(Station s : stationList.values()) {
        	s.setUnmarked();
            if(s.equals(stationList.get(origin))) {
                distance.put(s,0);
                foundStations.put(s,0);
                if(!getMinDistFromOrigin.containsKey(origin))
                getMinDistFromOrigin.put(origin,new TreeMap<String,Integer>());
               
            }
            else {
                distance.put(s,Integer.MAX_VALUE);
            }
            q.add(s);
        }
        
        //Algorithm
        while(!q.isEmpty()) {
            Station u = nearestStation(foundStations);
            
            q.remove(u);
            foundStations.remove(u);
            u.setMarked();

            Iterator<Station> iter = u.getAdjacentStations().keySet().iterator();
            while(iter.hasNext()) {
                Station neighbour= iter.next();
                if(q.contains(neighbour)) {
                    int alt = distance.get(u) + u.getAdjacentStations().get(neighbour);
                    if(alt < distance.get(neighbour)) {
                        distance.put(neighbour,alt);
                        foundStations.put(neighbour,alt);
                        
                        //Put Distance source->Dest in the lookup table
                        getMinDistFromOrigin.get(origin).put(neighbour.getName(),distance.get(neighbour));
                        //Put distance Dest-->Source in the lookup table
                        if(!getMinDistFromOrigin.containsKey(neighbour.getName())) getMinDistFromOrigin.put(neighbour.getName(),new TreeMap<String,Integer>());
                        getMinDistFromOrigin.get(neighbour.getName()).put(origin,distance.get(neighbour));
                       
                    }            
                }
            }
        }
        
        
}
   
   
   
   
    /**
     * Given a route between two stations, compute the total distance
     * of this route.
     *
     * @param path  the list of stations in the route (as String objects)
     * @return      the length of the route between the first station
     *              and last station in the list   
     */
 
    public int optimisedDistance(String origin,String destination) {
    	if(!getMinDistFromOrigin.containsKey(origin)|| !getMinDistFromOrigin.get(origin).containsKey(destination)) {
    		dijkstra(origin);
    	}

    	return getMinDistFromOrigin.get(origin).get(destination);
        }
	
	/**
	 * Given a route between two stations, compute the total distance 
	 * of this route.
	 * 
	 * @param path	the list of stations in the route (as String objects)
	 * @return		the length of the route between the first station
	 * 				and last station in the list	
	 */
	public int findTotalDistance(ArrayList<String> path)
	{
		if (path == null || path.isEmpty())
			return 0;

		int distance = 0;

		for (int stationNo = 0; stationNo < path.size() - 1; stationNo++)
		{
			Station current = stationList.get(path.get(stationNo));
			Station next = stationList.get(path.get(stationNo + 1));

			distance += current.getAdjacentStations().get(next);
		}

		return distance;
	}
	

	
	/**
	 * Return the ratio between the length of the shortest route between two
	 * stations (in terms of distance) and the actual distance between the 
	 * two stations (computed using computeDistance())
	 * 
	 * In other words, 
	 * let d1 = distance of shortest route between the two stations as computed
	 *          by the method routeMinStop() (from Stage 1).
	 * let d2 = distance between two stations as computed by the method
	 *          computeDistance() 
	 *          
	 * The method returns d1/d2 (as a double)
	 * 
	 * @param origin		the starting station
	 * @param destination 	the ending station
	 * @return	s			the ratio d1/d2 as explained above
	 */
	 public double computeRatio(String origin, String destination) {
	        return (double) optimisedDistance(origin,destination)/computeDistance(origin,destination);
	    }
	   
	
	/**
	 * Return the ratio as computed by computeRatio() method for all 
	 * pairs of station in the rail network.
	 * 
	 * The ratios should be stored in a HashMap<String,HashMap<String,Double>>,
	 * that is, the ratio between station a and b can be obtained by calling
	 * 
	 *    computeAllRatio().get(a).get(b)
	 * 
	 * @return a hashmap containing the ratios
	 */
	public HashMap<String,HashMap<String,Double>> computeAllRatio() {
		HashMap<String,HashMap<String,Double>> hMap = new HashMap<>();
        ArrayList<String> k = new ArrayList<String>(stationList.keySet()); // get all stations in a list
        for(int i=0;i<k.size();i++) {
        	
            HashMap<String,Double> h = null;
            if(!hMap.containsKey(k.get(i))) h= new HashMap<String,Double>();
            else h = hMap.get(k.get(i));
            for(int j=i+1;j<k.size();j++) {
            	HashMap<String,Double> z = null;
            	if(!hMap.containsKey(k.get(j))) z = new HashMap<String,Double>();
                else z = hMap.get(k.get(j));
                
            	//compute ratio of between station i and j
            	double d = computeRatio(k.get(i),k.get(j));
                h.put(k.get(j),d);
                z.put(k.get(i),d);
                //ratio (i,j) == ratio(j,i)
                hMap.put(k.get(i),h);
                hMap.put(k.get(j),z);
            }
        }
        return hMap;}
	
	/**
	 * The method finds the shortest route (in terms of number of stops)
	 * between the origin station and the destination station, taking 
	 * into account the available routes in the rail network.
	 * 
	 * The route is returned as an ArrayList<String> containing the lines taken,
	 * any transfer between lines, and the names of the stations on each line,
	 * including the origin and the destination stations.
	 * 
	 * Please see the assignment specification for more information.
	 * 
	 * If the route cannot be completed (there is no path between origin and
	 * destination), then return an empty ArrayList<String>
	 * 
	 * If the origin or the destination stations are not in the list of stations,
	 * return an empty ArrayList<String>. 
	 * 
	 * If the origin and the destination stations are the same, return an 
	 * ArrayList<String> containing the station.
	 * 
	 * @param origin		the starting station
	 * @param destination	the end station
	 * @return				the route taken
	 */
	public ArrayList<String> routeMinStopWithRoutes(String origin, String destination)
    {
		if (!stationList.containsKey(origin) || !stationList.containsKey(destination))
			return new ArrayList<>();

		//clear data sets
		pathFromOrigin.clear();
		distanceFromOrigin.clear();

		Station originStation = stationList.get(origin);
		Station destinationStation = stationList.get(destination);
		String currentLine = ""; //initial line value is empty

		//add original station
		distanceFromOrigin.put(originStation, 0);

		//launch recursion, true if track been found
		if (routeMinStopWithRoutesRec(originStation, destinationStation, currentLine))
			return pathFromOrigin.get(destinationStation);

		//if recursion returned false, nothing was found, return empty list
		return new ArrayList<>();
	}

	private Boolean routeMinStopWithRoutesRec(Station origin, Station destination, String currentLine)
	{
		//if we're on dest already, return true
		if (origin.equals(destination))
			return true;

		origin.setMarked();

		//get list of train lines the origin station is in
		ArrayList<String> lines = origin.getTrainLines();

		//iterate through the train line list
		for(String line : lines)
		{
			//0 - one station forward, 1 - one station backwards
			Station[] neighbours = new Station[2];

			//put neighbour stations indexes, if out of bound will remain null
			if (origin.getNumberByLine(line) + 1 <= lineList.get(line).getStationCount())
				neighbours[0] = stationList.get(lineList.get(line).getStation(origin.getNumberByLine(line) + 1));

			if (origin.getNumberByLine(line) - 1 >= 1)
				neighbours[1] = stationList.get(lineList.get(line).getStation(origin.getNumberByLine(line) - 1));

			for (Station neighbour : neighbours)
			{
			    //if neighbour exists and if either not in the list or distance to it would be better
				if (neighbour != null && (!distanceFromOrigin.containsKey(neighbour) ||
                        distanceFromOrigin.get(neighbour) > distanceFromOrigin.get(origin) + 1))
				{
					neighbour.setUnmarked();

					//this condition will only apply for the initial call
					ArrayList<String> list;
					if (pathFromOrigin.containsKey(origin))
					{
						list = new ArrayList<>(pathFromOrigin.get(origin));
					}
					else
					{
						list = new ArrayList<>();
					}

					//if changing a train line, add direction to the list
					if (currentLine.compareToIgnoreCase(line) != 0)
					{
						list.add(getLineDirection(line, origin, neighbour));
						list.add(origin.getName());
					}
					list.add(neighbour.getName());  //add this neighbour to the list

                    //fill both path and distance lists with updated/new data
					pathFromOrigin.put(neighbour, list);
					distanceFromOrigin.put(neighbour, distanceFromOrigin.get(origin) + 1);
				}
			}
		}

		Boolean result = false;

		//call recursively both neighbours for each train line station is in
		for(String line : lines)
		{
			//0 - one station forward, 1 - one station backwards
			Station[] neighbours = new Station[2];

			//same code as above, just get indexes of the stations
			if (origin.getNumberByLine(line) + 1 <= lineList.get(line).getStationCount())
				neighbours[0] = stationList.get(lineList.get(line).getStation(origin.getNumberByLine(line) + 1));

			if (origin.getNumberByLine(line) - 1 >= 1)
				neighbours[1] = stationList.get(lineList.get(line).getStation(origin.getNumberByLine(line) - 1));

			for (Station neighbour : neighbours)
			{
			    //call only if exists and not marked
				if (neighbour != null && !neighbour.isMarked())
					result |= routeMinStopWithRoutesRec(neighbour, destination, line); //if true came from any iteration, result would be true
			}
		}
		return result;
	}

	/**
     * Given two stations in a order (current, next) and a train line
     * these stations are on. Returns direction in a form
     * <Line Name> towards <Line End> from <Line Start>
	 *
	 * @param line code of line the stations are on
	 * @param origin current station
	 * @param next next station
	 * @return string with correct direction
	 */
	private String getLineDirection(String line, Station origin, Station next)
	{
		int originIndex = origin.getNumberByLine(line);
		int nextIndex = next.getNumberByLine(line);

		TrainLine trainLine = lineList.get(line);
		String result;

		if (originIndex < nextIndex)
		{
			result = trainLine.getName() + " towards " + trainLine.getEnd() + " from " + trainLine.getStart();
		}
		else
		{
			result = trainLine.getName() + " towards " + trainLine.getStart() + " from " + trainLine.getEnd();
		}

		return result;
	}
	
	public static void main(String[] args)
	{
		String stationData = "src/data/station_data.csv";
		String connectionData = "src/data/adjacent_stations.csv";
		String linesData = "src/data/lines_data.csv";
		RailNetworkAdvanced r = new RailNetworkAdvanced(stationData,connectionData,linesData);

		System.out.println(r.routeMinStopWithRoutes( "Blacktown", "Parramatta"));
	}
}
