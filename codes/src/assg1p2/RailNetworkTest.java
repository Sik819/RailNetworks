package assg1p2;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.TreeSet;

import org.junit.Before;
import org.junit.Test;

public class RailNetworkTest {

    RailNetworkAdvanced r;
    @Before public void initialize() {
        String stationData = "src/data/station_data.csv";
        String connectionData = "src/data/adjacent_stations.csv";
        String lineData = "";
        r = new RailNetworkAdvanced(stationData,connectionData, lineData);
    }

    /** Tests for routeMinDistance with no failures **/
    @Test
    public void routeMinDistanceTest1() {
        String origin = "Hornsby";
        String destination = "Hornsby";
        String[] expected = {"Hornsby"};

        ArrayList<String> actual = r.routeMinDistance(origin, destination);
        assertArrayEquals(expected, actual.toArray());
    }


    @Test
    public void routeMinDistanceTest2() {
        String origin = "Hornsby";
        String destination = "Chatswood";
        String[] expected = {"Hornsby","Waitara","Wahroonga","Warrawee",
                "Turramurra","Pymble","Gordon","Killara",
                "Lindfield","Roseville","Chatswood"};

        ArrayList<String> actual = r.routeMinDistance(origin, destination);
        assertArrayEquals(expected, actual.toArray());
    }

    @Test
    public void routeMinDistanceTest3() {
        String origin = "Wolli Creek";
        String destination = "Town Hall";
        String[] expected = {"Wolli Creek","Tempe","Sydenham","St Peters",
                "Erskineville","Redfern","Central","Town Hall"};


        ArrayList<String> actual = r.routeMinDistance(origin, destination);
        assertArrayEquals(expected, actual.toArray());
    }

    @Test
    public void routeMinDistanceTest4() {
        String origin = "East Richmond";
        String destination = "Hurstville";
        String[] expected = {"East Richmond","Clarendon","Windsor","Mulgrave",
                "Vineyard","Riverstone","Schofields","Quakers Hill","Marayong",
                "Blacktown","Seven Hills","Toongabbie","Pendle Hill",
                "Wentworthville","Westmead","Parramatta","Harris Park",
                "Granville","Clyde","Auburn","Lidcombe","Flemington","Homebush",
                "Strathfield","Burwood","Croydon","Ashfield","Summer Hill",
                "Lewisham","Petersham","Stanmore","Newtown","Macdonaldtown",
                "Redfern","Erskineville","St Peters","Sydenham","Tempe",
                "Wolli Creek","Arncliffe","Banksia","Rockdale","Kogarah",
                "Carlton","Allawah","Hurstville"};


        ArrayList<String> actual = r.routeMinDistance(origin, destination);
        assertArrayEquals(expected, actual.toArray());
    }

    /** Test for routeMinStop with no failures **/
    @Test
    public void routeMinStopTest1() {
        String origin = "Hornsby";
        String destination = "Hornsby";
        String[] expected = {"Hornsby"};

        ArrayList<String> actual = r.routeMinDistance(origin, destination);
        assertArrayEquals(expected, actual.toArray());
    }

    @Test
    public void routeMinStopTest2() {
        String origin = "Hornsby";
        String destination = "Epping";
        String[] expected = {"Hornsby","Normanhurst","Thornleigh","Pennant Hills",
                "Beecroft","Cheltenham","Epping"};

        ArrayList<String> actual = r.routeMinStop(origin, destination);
        assertArrayEquals(expected, actual.toArray());
    }

    @Test
    public void routeMinStopTest3() {
        String origin = "Wolli Creek";
        String destination = "Town Hall";
        String[] expected = {"Wolli Creek","International Airport",
                "Domestic Airport","Mascot","Green Square",
                "Central","Town Hall"};


        ArrayList<String> actual = r.routeMinStop(origin, destination);
        assertArrayEquals(expected, actual.toArray());
    }

    @Test
    public void routeMinStopTest4() {
        String origin = "Richmond";
        String destination = "Hurstville";
        String[] expected = {"Richmond","East Richmond","Clarendon","Windsor","Mulgrave",
                "Vineyard","Riverstone","Schofields","Quakers Hill","Marayong",
                "Blacktown","Seven Hills","Toongabbie","Pendle Hill",
                "Wentworthville","Westmead","Parramatta","Harris Park",
                "Granville","Clyde","Auburn","Lidcombe","Berala","Regents Park",
                "Birrong","Yagoona","Bankstown","Punchbowl","Wiley Park",
                "Lakemba","Belmore","Campsie","Canterbury","Hurlstone Park",
                "Dulwich Hill","Marrickville","Sydenham","Tempe","Wolli Creek",
                "Arncliffe","Banksia","Rockdale","Kogarah","Carlton","Allawah",
                "Hurstville"};

        ArrayList<String> actual = r.routeMinStop(origin, destination);
        assertArrayEquals(expected, actual.toArray());
    }

    /** Tests for findTotalDistance **/
    @Test
    public void findTotalDistanceTest1() {
        String origin = "Hornsby";
        String destination = "Epping";

        int expected = 9703;
        int actual = r.findTotalDistance(r.routeMinDistance(origin, destination));
        assertEquals(expected,actual);
    }

    @Test
    public void findTotalDistanceTest2() {
        String origin = "East Richmond";
        String destination = "Hurstville";

        int expected = 70768;
        int actual = r.findTotalDistance(r.routeMinDistance(origin, destination));
        assertEquals(expected,actual);
    }
}