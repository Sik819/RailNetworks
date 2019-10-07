package assg1p2;

public class TrainLine {

    private String code;
    private String name;
    private String start;
    private String end;
    private int stationCount;
    private String[] stations;

    public TrainLine(String code, String name, String start, String end, int stationCount)
    {
        this.code = code;
        this.name = name;
        this.start = start;
        this.end = end;
        this.stationCount = stationCount;
        this.stations = new String[stationCount + 1];
    }

    /* Getters */
    public String getCode() { return this.code; }
    public String getName() { return this.name; }
    public String getStart() { return this.start; }
    public String getEnd() { return this.end; }
    public int getStationCount() { return this.stationCount; }
    public String[] getStations() { return this.stations; }
    public String getStation(int index) { return this.stations[index]; }


    /*  */
    public void addStation(int index, String station)
    {
        //System.out.println(this.code + " - " + station + " " + index);
        this.stations[index] = station;
    }

    /* compareTo method since we are implementing Comparable */
    public int compareTo(TrainLine o) {
        if (o != null)
            return this.name.compareTo(o.getName());
        else
            return 1;
    }

    /*
     * Prints out the station name, followed by its neighbours and distances
     * to them
     */
    @Override
    public String toString()
    {
        return this.code + " {from " + this.start + " to " + this.end + "}";
    }
}
