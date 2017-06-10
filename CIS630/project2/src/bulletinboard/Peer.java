package bulletinboard;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.net.DatagramPacket;
import java.net.DatagramSocket;
import java.net.InetAddress;
import java.net.SocketException;
import java.net.SocketTimeoutException;
import java.net.UnknownHostException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Random;
import java.util.TimeZone;
import java.util.TreeMap;

/**
 * Class representing a single Peer in the peer-to-peer ring based system.
 *
 * @author abhishek
 */
public class Peer {
	private static final String PROBE = "PROBE";
	private static final String PROBE_REPLY = "PROBE-REPLY";
	private static final String ELECTION = "ELECTION";
	private static final String ELECTED = "ELECTED";
	private static final String TOKEN = "TOKEN";
	private static final long TIMEOUT = 1000;
	private static final long SLEEP = 1000;

	private static SimpleDateFormat formatter = new SimpleDateFormat("mm:ss");
	static {
		formatter.setTimeZone(TimeZone.getTimeZone("UTC"));
	}

	private static final String INPUT = "-i";
	private static final String CONFIG = "-c";
	private DatagramSocket socket;
	private int port;
	private int rangeStart, rangeEnd;
	private int nextHop, lastHop;
	private int stateCount = 0;
	private long goLiveTime;
	private long join, leave;
	private boolean alive;
	private boolean probeSuccess;
	private Receive receiveThread;
	private Probe probingThread;
	private NavigableMap<Long, String> messages;

	public long getJoin() {
		return join;
	}

	private String getTimestamp() {
		return "[" + formatter.format(new Date(System.currentTimeMillis() - goLiveTime)) + "]";
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

	private void reInitPorts() {
		this.nextHop = -1;
		this.lastHop = -1;
	}

	public Peer(int port, int rangeStart, int rangeEnd, long join, long leave) {
		this.port = port;
		this.rangeStart = rangeStart;
		this.rangeEnd = rangeEnd;
		this.join = join;
		this.leave = leave;
		this.nextHop = -1;
		this.lastHop = -1;
		System.out.println("Port : " + this.port);
		this.setAlive(false);
		this.probeSuccess = false;
	}

	public void send(String message, int destPort) throws PeerException {
		if (!isAlive())
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

	class Probe extends Thread {

		public Probe() {
			reInitPorts();
		}

		@Override
		public void run() {
			System.out.println("Started Probing " + stateCount);
			probeSuccess = false;
			for (int probePort = port + 1; probePort <= rangeEnd; probePort++) {
				if (probePort == port) {
					if (probePort == rangeEnd)
						probePort = rangeStart - 1;
					continue;
				}
				try {
					send(PROBE, probePort);
					Thread.sleep(SLEEP / 10);
				} catch (PeerException | InterruptedException e) {
					System.err.println("Error occured while sending probing message");
					System.exit(1);
				}
				if (probeSuccess) {
					nextHop = probePort;
					break;
				}
				if (probePort == rangeEnd)
					probePort = rangeStart - 1;
			}
			/*
			 * Initiate probing after you have received a successful probe
			 */
			if (probeSuccess) {
				try {
					sleep(new Random().nextInt(1000));
				} catch (InterruptedException e) {
					System.err.println("Error occurred while initializing election");
					System.exit(1);
				}
				try {
					send(ELECTION + ":" + port, nextHop);
				} catch (PeerException e) {
					System.err.println("Error occurred while sending election");
					System.exit(1);
				}
			}
		}
	};

	class Receive extends Thread {
		@Override
		public void run() {
			System.out.println("Receive listener is up");
			while (Peer.this.isAlive()) {
				try {
					byte buffer[] = new byte[256];
					DatagramPacket packet = new DatagramPacket(buffer, buffer.length);
					socket.receive(packet);
					int senderPort = packet.getPort();
					String message = new String(packet.getData(), 0, packet.getLength());
					System.out.println(getTimestamp() + " : " + message + " sender port : " + senderPort
							+ " next hop : " + nextHop + " last hop :" + lastHop);
					if (message.equals(PROBE)) {
						if (lastHop == -1 || (Math.floorMod(senderPort - port, rangeEnd - rangeStart) > Math
								.floorMod(lastHop - port, rangeEnd - rangeStart))) {
							lastHop = senderPort;
							System.out.println(getTimestamp() + " : " + "Changed previous hop " + lastHop);
							send(PROBE_REPLY, senderPort);
						}
					} else if (message.equals(PROBE_REPLY)) {
						if (nextHop == -1 || (Math.floorMod(senderPort - port, rangeEnd - rangeStart) < Math
								.floorMod(nextHop - port, rangeEnd - rangeStart))) {
							nextHop = senderPort;
							probeSuccess = true;
							System.out.println(getTimestamp() + " : " + "Changed next hop " + nextHop);
						}
					}
					/*
					 * If probing was not successful and last hope was not
					 * initialized / or last hop was not same as the sender,
					 * continue
					 */
					if (!probeSuccess || senderPort != lastHop)
						continue;

					if (message.startsWith(ELECTION)) {

						int electionPort = Integer.parseInt((message.split(":"))[1]);
						if (electionPort < port) {
							// Replace the
							send(ELECTION + ":" + port, nextHop);
						} else if (electionPort > port) {
							// Forward the election message as is
							send(message, nextHop);
						} else {
							send(ELECTED + ":" + port, nextHop);
						}

					} else if (message.startsWith(ELECTED)) {
						int electionPort = Integer.parseInt((message.split(":"))[1]);
						if (electionPort < port) {

						} else if (electionPort > port) {
							send(message, nextHop);
						} else {
							generateToken();
							sendMessages();
						}
					} else if (message.startsWith(TOKEN)) {
						send(message, nextHop);
					}
				} catch (SocketTimeoutException e) {
					if (!probingThread.isAlive()) {
						probingThread = new Probe();
						probingThread.start();
					}
				} catch (PeerException | IOException e) {
					System.err.println("There was a problem while " + "receiving server request :" + e.getMessage()
							+ " " + e.getClass().getName());
					socket.close();
					System.exit(1);
				}
			}
		}
	};

	public void receive() {
		receiveThread = new Receive();
		receiveThread.start();
	};

	public void sendMessages() throws PeerException {
		send(TOKEN + ":" + port, nextHop);
	}

	public int generateToken() {
		return new Random().nextInt(10000);
	}

	public void probe() {
		probingThread = new Probe();
		probingThread.start();
	};

	public boolean isAlive() {
		return alive;
	}

	public void setAlive(boolean alive) {
		this.alive = alive;
	}

	@SuppressWarnings("deprecation")
	private static void schedule(Peer peer) throws PeerException {
		System.out.println("Scheduling peer : " + peer.port);
		long goLive = peer.getGoLiveTime();
		long join = peer.getJoin();
		long leave = peer.getLeave();
		System.out.println("waiting to join");
		while ((System.currentTimeMillis()) - goLive < join) {
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
		while ((System.currentTimeMillis()) - goLive < leave) {
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
		if (reader != null)
			buffReader = new BufferedReader(reader);
		else {
			System.err.println("Error occured while reading config file");
			System.exit(1);
		}
		String line;
		while ((line = buffReader.readLine()) != null) {
			if (line.trim().isEmpty())
				continue;
			String[] split = line.split("\t");
			long timestamp = 0;
			timestamp = formatter.parse(split[0].trim()).getTime();
			messages.put(timestamp, split[1].trim());
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
		if (reader != null)
			buffReader = new BufferedReader(reader);
		else {
			System.err.println("Error occured while reading config file");
			System.exit(1);
		}
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
		if (peer == null || messages == null) {
			System.err.println("System might be in an inconsistent state. "
					+ "Either the peer or the messages were not initialized");
			System.exit(1);

		}
		peer.setGoLiveTime(System.currentTimeMillis());
		peer.setMessages(messages);
		try {
			schedule(peer);
		} catch (PeerException e) {
			System.err.println("Error occurred while scheduling peer.");
			System.exit(1);
		}
	}
};

class PeerException extends Exception {

	private static final long serialVersionUID = -1165295127240695115L;

	public PeerException(String message) {
		super(message);
	}
};
