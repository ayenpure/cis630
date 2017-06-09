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
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Random;
import java.util.TimeZone;
import java.util.TreeMap;

/*
 * TODO : If you are the elected, when you receive an election message, reply
 * if you receive an elected message, just proceed sending messages.
 * Make sure the message you had already sent is received before sending 
 * out a new message.
 * If I have the token right now and get the election message with a higher port,
 * then quit sending messages, give up token, and forward election message.
 * else drop forwarding the election message
 */

/**
 * Class representing a single Peer in the peer-to-peer ring based system.
 * 
 * @author abhishek
 */
public class Peer {
	private static final String TOKEN = "TOKEN";
	private static final String PROBE = "PROBE";
	private static final String ELECTION = "ELECTION";
	private static final String ELECTED = "ELECTED";
	private static final String PROBE_REPLY = "PROBE-REPLY";
	private static final String MESSAGE = "MESSAGE";
	private static final long TIMEOUT = 1000;
	private static final long SLEEP = 1000;
	private static SimpleDateFormat formatter = new SimpleDateFormat("mm:ss");
	static {
		formatter.setTimeZone(TimeZone.getTimeZone("UTC"));
	}
	
	private String getTimeStamp() {
		return "[" + formatter.format(new Date(System.currentTimeMillis() - goLiveTime)) + "] : ";
	}

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
	private boolean alive;
	private boolean probeSuccess;
	private boolean hasToken;
	private boolean tokenChanged;
	private boolean forwarding;
	private int currentToken;
	private long tokenTimeStamp;
	private Receive receiveThread;
	private Probe probingThread;
	private NavigableMap<Long, String> messages;

	private void reInitToken() {
		hasToken = false;
		tokenChanged = true;
		currentToken = -1;
		tokenTimeStamp = -1L;
	}

	private void reInitToken(String token, long receivedTimeStamp) {
		int previousToken = currentToken;
		hasToken = true;
		currentToken = Integer.parseInt((token.split(":"))[1]);
		tokenTimeStamp = receivedTimeStamp;
		if(previousToken != currentToken)
			System.out.println(getTimeStamp()+"token ["+currentToken+"] was received");
	}

	private void generateToken(long receivedTimestamp) {
		hasToken = true;
		tokenChanged = true;
		tokenTimeStamp = receivedTimestamp;
		currentToken = Math.abs(new Random(receivedTimestamp).nextInt() % 10000);
		System.out.println(getTimeStamp() +"new token generated ["+currentToken+"]");
	}

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

	private void reInitPorts() {
		this.nextHop = Integer.MAX_VALUE;
		this.previousHop = -1;
		forwarding = false;
	}

	public Peer(int port, int rangeStart, int rangeEnd, long join, long leave) {
		this.port = port;
		this.rangeStart = rangeStart;
		this.rangeEnd = rangeEnd;
		this.join = join;
		this.leave = leave;
		this.hasToken = false;
		this.currentToken = -1;
		this.forwarding = false;
		this.setAlive(false);
		this.probeSuccess = false;
		this.reInitPorts();
		this.reInitToken();
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

		@Override
		public void run() {
			probeSuccess = false;
			for (int i = port + 1; i <= rangeEnd && probeSuccess == false; i++) {
				if (i == port) {
					if (i == rangeEnd)
						i = rangeStart - 1;
					continue;
				}
				try {
					send(PROBE, i);
					Thread.sleep(SLEEP);
					if (probeSuccess) {
						nextHop = i;
						Thread.sleep(SLEEP / 2);
						System.out.println(getTimeStamp() + "started election, send election message to client ["+nextHop+"]");
						send(ELECTION + ":" + port, nextHop);
						break;
					}
				} catch (PeerException | InterruptedException e) {
					System.exit(1);
				}
				if (i == rangeEnd)
					i = rangeStart - 1;
			}
		}
	};

	private void sendMessageFromInput(long tokenTimeStamp) throws PeerException {
		/*
		 * if there is a message, forward it else forward token
		 */
		NavigableMap<Long, String> headMap = messages.headMap(tokenTimeStamp, true);
		if (headMap.isEmpty()) {
			send(TOKEN + ":" + currentToken, nextHop);
			if(tokenChanged) {
				tokenChanged = false;
				System.out.println(getTimeStamp()+"token ["+currentToken+"] was sent to client ["+nextHop+"]");	
			}
		}
		else {
			Long firstKey = headMap.firstKey();
			String message = headMap.get(firstKey);
			messages.remove(firstKey);
			send(MESSAGE + ":" + currentToken + ":" + message, nextHop);
		}
	}

	class Receive extends Thread {

		@Override
		public void run() {
			while (Peer.this.isAlive()) {
				try {
					byte buffer[] = new byte[256];
					DatagramPacket packet = new DatagramPacket(buffer, buffer.length);
					socket.receive(packet);
					long receivedTimestamp = System.currentTimeMillis();
					int senderPort = packet.getPort();
					String message = new String(packet.getData(), 0, packet.getLength());

					// Probing routine
					if (message.equals(PROBE)) {
						// Check for cyclic greater
						if (senderPort > previousHop) {
							previousHop = senderPort;
							System.out.println(getTimeStamp() + "previous hop is changed to client [" + previousHop + "]");
							send(PROBE_REPLY, senderPort);
						}
					} else if (message.equals(PROBE_REPLY)) {
						// Check for cyclic smaller
						if (senderPort < nextHop) {
							nextHop = senderPort;
							System.out.println(getTimeStamp() + "next hop is changed to client [" + nextHop + "]");
							probeSuccess = true;
						}
					}

					if (senderPort != previousHop || nextHop == Integer.MAX_VALUE)
						continue;
					// Election routine
					if (message.startsWith(ELECTION)) {
						forwarding = false;
						String[] split = message.split(":");
						int electionPort = Integer.parseInt(split[1].trim());
						if (electionPort < port) {
							send(ELECTION + ":" + port, nextHop);
							System.out.println(getTimeStamp() + "relayed election message, replaced leader");
						} else if (electionPort == port) {
							send(ELECTED + ":" + port, nextHop);
						} else {
							reInitToken();
							send(message, nextHop);
							System.out.println(getTimeStamp() + "relayed election message, leader: client ["+electionPort+"]");
						}
					} else if (message.startsWith(ELECTED)) {
						forwarding = true;
						String[] split = message.split(":");
						int electedPort = Integer.parseInt(split[1].trim());
						if (electedPort == port) {
							/*
							 * Generate token Send messages until now Forward
							 * Token
							 */
							if(!hasToken) {
								System.out.println(getTimeStamp() + "leader selected");
								generateToken(receivedTimestamp);
								sendMessageFromInput(tokenTimeStamp);
							} else {
								sendMessageFromInput(tokenTimeStamp);
							}
						} else if (electedPort > port) {
							/*
							 * Forward the elected message until it is received
							 * by the owner and the current port is not the owner
							 */
							send(message, nextHop);
						} else {
							sendMessageFromInput(tokenTimeStamp);
						}
					}

					if (!forwarding)
						continue;
					// Message passing after receiving token
					if (message.startsWith(TOKEN)) {
						reInitToken(message, receivedTimestamp);
						/*
						 * Send messages until now Forward token
						 */
						sendMessageFromInput(tokenTimeStamp);
					} else if (message.startsWith(MESSAGE)) {
						if (hasToken) {
							/*
							 * Check if the received message belongs to the
							 * token with us.
							 */
							String[] split = message.split(":", 3);
							if (Integer.parseInt(split[1].trim()) == currentToken) {
								sendMessageFromInput(receivedTimestamp);
							}
						} else
							send(message, nextHop);
					}
				} catch (SocketTimeoutException e) {
					reInitPorts();
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
		long goLive = peer.getGoLiveTime();
		long join = peer.getJoin();
		long leave = peer.getLeave();
		while ((System.currentTimeMillis()) - goLive < join) {
			try {
				Thread.sleep(SLEEP);
			} catch (InterruptedException e) {
				System.err.println("Sleep Error");
				System.exit(1);
			}
		}
		System.out.println(peer.getTimeStamp() + "joining");
		peer.setAlive(true);
		peer.initSocket();
		peer.receive();
		peer.probe();
		while ((System.currentTimeMillis()) - goLive < leave) {
			try {
				Thread.sleep(SLEEP);
			} catch (InterruptedException e) {
				System.err.println("Sleep Error");
				System.exit(1);
			}
		}
		System.out.println(peer.getTimeStamp() + "leaving");
		peer.setAlive(false);
		peer.receiveThread.stop();
		peer.probingThread.stop();
		peer.socket.close();
	}

	private void initSocket() throws PeerException {
		try {
			socket = new DatagramSocket(this.port);
			socket.setSoTimeout(2 * (int) TIMEOUT);
		} catch (SocketException e) {
			throw new PeerException("Failed to create socket with port " + port);
		}

	}

	private static NavigableMap<Long, String> getMessageMap(String fileName) throws IOException, ParseException {
		NavigableMap<Long, String> messages = new TreeMap<Long, String>();
		/*File input;
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
			long timestamp = formatter.parse(split[0].trim()).getTime();
			messages.put(timestamp, split[1].trim());
		}*/
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
			System.err.println("System might be in an inconsistent state. "
					+ "Either the peer or the messages were not initialized");
			System.exit(1);
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