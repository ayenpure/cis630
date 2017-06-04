package bullitenboard;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.net.DatagramPacket;
import java.net.DatagramSocket;
import java.net.InetAddress;
import java.net.SocketException;
import java.net.UnknownHostException;
import java.util.HashMap;
import java.util.Map;

import javax.swing.text.html.HTMLDocument.HTMLReader.SpecialAction;

/**
 * Class representing a single Peer in the peer-to-peer ring based system.
 * 
 * @author abhishek
 */
public class Peer {
	// private static final String TOKEN = "TOKEN";

	private static final String CONFIG = "-c";
	private DatagramSocket socket;
	private int rangeStart, rangeEnd;
	private int nextHop, previousHop;
	private int port;
	private boolean alive;

	public Peer(int port, int rangeStart, int rangeEnd) throws PeerException {
		this.port = port;
		this.rangeStart = rangeStart;
		this.rangeEnd = rangeEnd;
		this.nextHop = (port == rangeStart) ? rangeEnd : rangeStart;
		this.previousHop = -1;
		System.out.println("Port : " + port);
		System.out.println("Next hop : " + nextHop);
		this.setAlive(false);
		try {
			socket = new DatagramSocket(port);
		} catch (SocketException e) {
			throw new PeerException("Failed to create socket with port " + port);
		}
	}

	class Receive extends Thread {
		@Override
		public void run() {
			System.out.println("Receive listener is up");
			while (Peer.this.isAlive()) {
				try {
					byte buffer[] = new byte[256];
					DatagramPacket packet = new DatagramPacket(buffer, buffer.length);
					socket.receive(packet);
					System.out
							.println("Received message reads : " + new String(packet.getData(), 0, packet.getLength()));
				} catch (IOException e) {
					System.err.println("There was a problem while receiving server request :" + e.getMessage());
					socket.close();
					System.exit(1);
				}
			}
			socket.close();
		}
	};

	public void send(String message) throws PeerException {
		byte buffer[] = message.getBytes();
		InetAddress address = null;
		try {
			address = InetAddress.getLocalHost();
		} catch (UnknownHostException e) {
			throw new PeerException("Failed to get address for the serving peer : " + e.getMessage());
		}
		DatagramPacket packet = new DatagramPacket(buffer, buffer.length, address, nextHop);
		try {
			socket.send(packet);
		} catch (IOException e) {
			throw new PeerException("Failed to send request to the serving peer" + e.getMessage());
		}
	};

	public void receive() {
		Receive receive = new Receive();
		receive.start();
	};

	public boolean isAlive() {
		return alive;
	}

	public void setAlive(boolean alive) {
		this.alive = alive;
	}

	public static void main(String[] args) {

		Map<String, String> argMap = new HashMap<String, String>();
		argMap.put(args[0], args[1]);
		argMap.put(args[2], args[3]);
		argMap.put(args[4], args[5]);

		Peer peer = null;
		File config;
		FileReader reader = null;
		BufferedReader buffReader;

		config = new File(argMap.get(CONFIG));
		try {
			reader = new FileReader(config);
		} catch (FileNotFoundException e) {
			System.err.println("Failed to read config file");
		}
		buffReader = new BufferedReader(reader);

		String line;
		try {
			int port = 0, rangeStart = 0, rangeEnd = 0;
			while ((line = buffReader.readLine()) != null) {
				if (line.contains("client_port")) {
					String[] split = line.split(":");
					String[] range = split[1].split("-");
					rangeStart = Integer.parseInt(range[0].trim());
					rangeEnd = Integer.parseInt(range[1].trim());
				} else if (line.contains("my_port")) {
					String[] split = line.split(":");
					port = Integer.parseInt(split[1].trim());
				} else if (line.contains("join_time")) {

				} else if (line.contains("leave_time")) {

				}
			}
			peer = new Peer(port, rangeStart, rangeEnd);
		} catch (NumberFormatException | IOException | PeerException e) {
			System.err.println("Error occured while reading config file" + e.getMessage());
			if (peer != null)
				peer.socket.close();
			System.exit(1);
		}
		try {
			Thread.sleep(10000);
		} catch (InterruptedException e) {
			System.err.println("Error occurred while waiting");
			peer.socket.close();
			System.exit(1);
		}
		peer.setAlive(true);
		peer.receive();
		try {
			peer.send("");
		} catch (PeerException e) {
			System.err.println("Error occured while reading config file" + e.getMessage());
			peer.socket.close();
			System.exit(1);
		}
	}
};

class PeerException extends Exception {
	public PeerException(String message) {
		super(message);
	}
};