package bulletinboard;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.net.DatagramPacket;
import java.net.DatagramSocket;
import java.net.InetAddress;
import java.net.SocketException;
import java.net.SocketTimeoutException;
import java.net.UnknownHostException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.HashMap;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TimeZone;
import java.util.TreeMap;

/**
 * Class representing a single Peer in the peer-to-peer ring based system.
 * 
 * @author abhishek
 */
public class Peer {
	private static final String TOKEN = "TOKEN";
	private static final String PROBE = "PROBE";
	private static final String PROBE_REPLY = "PROBE-REPLY";
	private static final long TIMEOUT = 1000;
	private static final long SLEEP = 1000;
	
	class Message implements Serializable {
		private static final long serialVersionUID = 2111001394711772492L;
		String type;
		Object content;
	}

	private static final String INPUT = "-i";
	private static final String CONFIG = "-c";
	private DatagramSocket socket;
	private int port;
	private int rangeStart, rangeEnd;
	private int nextHop, previousHop;
	private long goLiveTime;
	private long join, leave;
	private long probetime;
	private boolean alive;
	private boolean probeSuccess;
	private Receive receiveThread;
	private Probe probingThread;
	private NavigableMap<Long, String> messages;

	public long getJoin() {
		return join;
	}
	
	public long getLeave() {
		return leave;
	}

	public long getGoLiveTime() {
		return goLiveTime;
	}

	public void setGoLiveTime(long goLiveTime) {
		
		this.goLiveTime = goLiveTime;
	}
	
	private void setMessages(NavigableMap<Long, String> messages) {
		this.messages = messages;
	}

	public Peer(int port, int rangeStart, int rangeEnd, long join, long leave) {
		this.port = port;
		this.rangeStart = rangeStart;
		this.rangeEnd = rangeEnd;
		this.join = join;
		this.leave = leave;
		this.nextHop = -1;
		this.previousHop = Integer.MAX_VALUE;
		System.out.println("Port : " + this.port);
		System.out.println("Next hop : " + nextHop);
		this.setAlive(false);
		this.probeSuccess = false;
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
					long receivedTimestamp = System.currentTimeMillis();
					int senderPort = packet.getPort();
					String message = new String(packet.getData(), 0, packet.getLength());
					//Check for cyclic great
					if(message.equals(PROBE)){
						if(senderPort < previousHop) {
							previousHop = senderPort;
							System.out.println("Previous Hop ===> " + previousHop);
						}
						send(PROBE_REPLY, senderPort);
					} else if (message.equals(PROBE_REPLY)) {
						//long elapsed = receivedTimestamp - probetime;
						if(senderPort > nextHop) {
							nextHop = senderPort;
							System.out.println("Next Hop ===> " + nextHop);
							probeSuccess = true;
						}
						
					}
				} catch(SocketTimeoutException e) {
						nextHop = -1;
						previousHop = -1;
						if(!probingThread.isAlive()) {
							probingThread = new Probe();
							probingThread.start();
						}
				} catch (PeerException | IOException e) {
					System.err.println("There was a problem while receiving server request :" + e.getMessage()
					+ " " + e.getClass().getName());
					socket.close();
					System.exit(1);
				}
			}
		}
	};
	
	class Probe extends Thread {
		@Override
		public void run() {
			System.out.println("Started Probing");
			probeSuccess = false;
			for(int i = port + 1; i <= rangeEnd; i++) {
				if(i == port) {
					if(i == rangeEnd)
						i = rangeStart - 1;
					continue;
				}
				//System.out.println("Probing " + i);
				try {
					send(PROBE, i);
					probetime = System.currentTimeMillis();
					Thread.sleep(SLEEP);
				} catch (PeerException | InterruptedException e) {
					//System.err.println("Error occured while probing client");
					System.exit(1);
				}
				if(probeSuccess) {
					//System.out.println("Probe successful with " + i);
					nextHop = i;
					break;
				}
				//System.out.println("Probe unsuccessful with " + i);
				if(i == rangeEnd)
					i = rangeStart - 1;
			}
		}
	};

	public void send(String message, int destPort) throws PeerException {
		if(!isAlive())
			return;
		byte buffer[] = message.getBytes();
		InetAddress address = null;
		try {
			address = InetAddress.getLocalHost();
		} catch (UnknownHostException e) {
			throw new PeerException("Failed to get address for the serving peer : " + e.getMessage());
		}
		DatagramPacket packet = new DatagramPacket(buffer, buffer.length, address, destPort);
		try {
			socket.send(packet);
		} catch (IOException e) {
			throw new PeerException("Failed to send request to the serving peer" + e.getMessage());
		}
	};

	public void receive() {
		receiveThread = new Receive();;
		receiveThread.start();
	};
	
	public void probe() {
		probingThread = new Probe();;
		probingThread.start();
	};

	public boolean isAlive() {
		return alive;
	}

	public void setAlive(boolean alive) {
		if(alive)
			System.out.println("Peer " + port + " is alive");
		this.alive = alive;
	}

	@SuppressWarnings("deprecation")
	private static void schedule(Peer peer) throws PeerException {
		System.out.println("Scheduling peer : " + peer.port);
		long goLive = peer.getGoLiveTime();
		long join = peer.getJoin();
		long leave = peer.getLeave();
		System.out.println("waiting to join");
		while((System.currentTimeMillis()) - goLive < join) {
			try {
				Thread.sleep(SLEEP);
			} catch (InterruptedException e) {
				System.err.println("Sleep Error");
				System.exit(1);
			}
			
		}
		System.out.println("joining");
		peer.setAlive(true);
		peer.initSocket();
		peer.receive();
		peer.probe();
		System.out.println("waiting to leave");
		while((System.currentTimeMillis()) - goLive < leave) {			
			try {
				Thread.sleep(SLEEP);
			} catch (InterruptedException e) {
				System.err.println("Sleep Error");
				System.exit(1);
			}
		}
		System.out.println("leaving");
		peer.setAlive(false);
		peer.receiveThread.stop();
		peer.probingThread.stop();
		peer.socket.close();
	}

	private void initSocket() throws PeerException {
		try {
			socket = new DatagramSocket(this.port);
			socket.setSoTimeout((int) TIMEOUT);
		} catch (SocketException e) {
			throw new PeerException("Failed to create socket with port " + port);
		}
		
	}

	private static NavigableMap<Long, String> getMessageMap(String fileName) throws IOException, ParseException {
		NavigableMap<Long, String> messages = new TreeMap<Long, String>();
		File input;
		FileReader reader = null;
		BufferedReader buffReader = null;
		input = new File(fileName);
		try {
			reader = new FileReader(input);
		} catch (FileNotFoundException e) {
			System.err.println("Failed to read config file");
			System.exit(1);
		}
		if(reader != null)
			buffReader = new BufferedReader(reader);
		else {
			System.err.println("Error occured while reading config file");
			System.exit(1);
		}
		SimpleDateFormat formatter = new SimpleDateFormat("mm:ss");
		formatter.setTimeZone(TimeZone.getTimeZone("UTC"));
		String line;
		while ((line = buffReader.readLine()) != null) {
			if(line.trim().isEmpty())
				continue;
			String[] split = line.split("\t");
			long timestamp = 0;
			timestamp = formatter.parse(split[0].trim()).getTime();
			messages.put(timestamp , split[1].trim());
		}	
		return messages;
	}

	private static Peer constructPeer(String fileName) {
		Peer peer = null;
		File config;
		FileReader reader = null;
		BufferedReader buffReader = null;

		config = new File(fileName);
		try {
			reader = new FileReader(config);
		} catch (FileNotFoundException e) {
			System.err.println("Failed to read config file");
			System.exit(1);
		}
		if(reader != null)
			buffReader = new BufferedReader(reader);
		else {
			System.err.println("Error occured while reading config file");
			System.exit(1);
		}
		SimpleDateFormat formatter = new SimpleDateFormat("mm:ss");
		formatter.setTimeZone(TimeZone.getTimeZone("UTC"));
		String line;
		try {
			int port = 0, rangeStart = 0, rangeEnd = 0;
			long join = 0, leave = 0;
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
					String[] split = line.split(":", 2);
					join = formatter.parse(split[1].trim()).getTime();
				} else if (line.contains("leave_time")) {
					String[] split = line.split(":", 2);
					leave = formatter.parse(split[1].trim()).getTime();
				}
			}
			peer = new Peer(port, rangeStart, rangeEnd, join, leave);
		} catch (NumberFormatException | IOException | ParseException e) {
			System.err.println("Error occured while reading config file" + e.getMessage());
			if (peer != null)
				peer.socket.close();
			System.exit(1);
		}
		return peer;
	}
		
	public static void main(String[] args) {

		Map<String, String> argMap = new HashMap<String, String>();
		argMap.put(args[0], args[1]);
		argMap.put(args[2], args[3]);
		argMap.put(args[4], args[5]);
		
		Peer peer = null;
		NavigableMap<Long, String> messages = null;
		try {
			peer = constructPeer(argMap.get(CONFIG));
			messages = getMessageMap(argMap.get(INPUT));
		} catch (IOException | ParseException e) {
			
		}
		if(peer == null || messages == null) {
			System.err.println("System might be in an inconsistent state. "
					+ "Either the peer of the messages were not initialized");
			System.exit(1);
		}
		peer.setGoLiveTime(System.currentTimeMillis());
		peer.setMessages(messages);
		try {
			schedule(peer);
		} catch (PeerException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

};

class PeerException extends Exception {
	
	private static final long serialVersionUID = -1165295127240695115L;

	public PeerException(String message) {
		super(message);
	}
};