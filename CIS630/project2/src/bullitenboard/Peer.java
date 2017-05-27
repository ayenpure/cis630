package bullitenboard;

import java.io.IOException;
import java.net.DatagramPacket;
import java.net.DatagramSocket;
import java.net.InetAddress;
import java.net.SocketException;
import java.net.UnknownHostException;

public class Peer {
	public static void main(String[] args) {
		DatagramSocket socket = null;
		try {
			socket = new DatagramSocket();
		} catch (SocketException e) {
			System.err.println("Failed to create requesting socket");
			System.exit(1);
		}
		byte[] buf = new byte[256];
        InetAddress address = null;
		try {
			address = InetAddress.getLocalHost();
		} catch (UnknownHostException e) {
			System.err.println("Failed to get address for the serving peer");
			System.exit(1);
		}
        DatagramPacket packet = new DatagramPacket(buf, buf.length, address, 4445);
        try {
			socket.send(packet);
		} catch (IOException e) {
			System.err.println("Failed to send request to the serving peer");
			System.exit(1);
		}
    
            // get response
        packet = new DatagramPacket(buf, buf.length);
        try {
			socket.receive(packet);
		} catch (IOException e) {
			System.err.println("Failed to receive response from the serving peer");
			System.exit(1);
		}

	    // display response
        String received = new String(packet.getData(), 0, packet.getLength());
        System.out.println("Quote of the Moment: " + received);
        socket.close();
	}
}
