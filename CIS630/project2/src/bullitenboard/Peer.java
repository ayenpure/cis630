package bullitenboard;

import java.io.IOException;
import java.net.DatagramPacket;
import java.net.DatagramSocket;
import java.net.InetAddress;
import java.net.SocketException;
import java.net.UnknownHostException;

/**
 * Class representing a single Peer in the peer-to-peer ring based system.
 * 
 * @author abhishek
 *
 */
public class Peer {

	private DatagramSocket socket;
	
	public Peer() {
		int sPort = 4445;
		try {
			socket = new DatagramSocket(sPort);
		} catch (SocketException e) {
			System.err.println("Failed to create socket with port " + sPort);
			System.exit(1);
		}
	}
	
	class Receive extends Thread {
		private boolean status = true;

		public void setStatus(boolean status) {
			this.status = status;
		}

		@Override
		public void run() {
			System.out.println("Client is up");
			byte buffer[] = new byte[256];
			while (status) {
				DatagramPacket packet = new DatagramPacket(buffer, buffer.length);
				try {
					socket.receive(packet);
					System.out.println("Client request reads : " + new String(packet.getData(), 0, packet.getLength()));
				} catch (IOException e) {
					System.err.println("Failed to receive request from client");
					System.exit(1);
				}
				String toSend = "Hello!! from the serving socket";
				InetAddress address = packet.getAddress();
				int dPort = packet.getPort();
				buffer = toSend.getBytes();
				packet = new DatagramPacket(buffer, buffer.length, address, dPort);
				try {
					socket.send(packet);
				} catch (IOException e) {
					System.err.println("Failed to respond to client request");
					System.exit(1);
				}
			}
		}
	};

	class Send extends Thread {
		@Override
		public void run() {
			String request = "Hello!! from the requesting socket";
			byte buffer[] = request.getBytes();
			InetAddress address = null;
			try {
				address = InetAddress.getLocalHost();
			} catch (UnknownHostException e) {
				System.err.println("Failed to get address for the serving peer");
				System.exit(1);
			}
			DatagramPacket packet = new DatagramPacket(buffer, buffer.length, address, 4445);
			try {
				socket.send(packet);
			} catch (IOException e) {
				System.err.println("Failed to send request to the serving peer");
				System.exit(1);
			}

			// get response
			packet = new DatagramPacket(buffer, buffer.length);
			try {
				socket.receive(packet);
			} catch (IOException e) {
				System.err.println("Failed to receive response from the serving peer");
				System.exit(1);
			}

			// display response
			String received = new String(packet.getData(), 0, packet.getLength());
			System.out.println("Quote of the Moment: " + received);
			//socket.close();
		};
	};

	public void send() {
		Send send = new Send();
		send.start();
	};
	
	public void receive() {
		Receive receive = new Receive();
		receive.start();
	};
	
	public static void main(String[] args) {
		Peer peer = new Peer();
		peer.receive();
		peer.send();
	};
};