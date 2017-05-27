package bullitenboard;

import java.io.IOException;
import java.net.DatagramPacket;
import java.net.DatagramSocket;
import java.net.InetAddress;
import java.net.SocketException;

/**
 * Class representing a single Peer in the peer-to-peer ring based system.
 * @author abhishek
 *
 */
public class ServingPeer extends Thread {

	private DatagramSocket socket;
	private boolean status = true;
	
	public void setStatus(boolean status) {
		this.status = status;
	}
	
	public ServingPeer() {
		int sPort = 4445;
		try {
			socket = new DatagramSocket(sPort);
		} catch (SocketException e) {
			System.err.println("Failed to create socket with port " + sPort);
			System.exit(1);
		}
	}
	
	
	
	@Override
	public void run() {
		System.out.println("Client is up");
		byte buffer[] = new byte[256];
		while(status){
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
	
	public static void main(String[] args) {
		ServingPeer peer = new ServingPeer();
		peer.start();
	}
}
