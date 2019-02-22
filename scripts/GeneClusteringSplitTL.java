import java.io.*;
import java.util.*;

public class GeneClusteringSplitTL {
  
  public static Hashtable htBBH = new Hashtable();
  public static Hashtable htSpe = new Hashtable();
  public static LinkedList lGene = new LinkedList();
  public static LinkedList lSpe = new LinkedList();
  public static String myspe = "";
  
  public static void main(String[] args) {
    String fileGene = args[1];
    String fileSpe = args[2];
    String dir = args[3];  
    String fileBBH = args[4];
    myspe = args[0];
    
    readGenes(fileGene);
    readSpe(fileSpe);
    readBBH(fileBBH);
    readGenomes(dir);
    
    int start = 0;
    int end = lGene.size();
    if (args.length>5)
      start = Integer.parseInt(args[5]);
    if (args.length>6)
      end = Integer.parseInt(args[6]);

System.out.println("Total number of genomes: "+(lSpe.size()-1)+" + 1");
//for (int i=0; i<lSpe.size(); i++)
//  System.out.print(lSpe.get(i)+" ");
//System.out.println();

    for (int i=start; i<end; i++) {
      String g1 = lGene.get(i).toString();
      Hashtable h1 = (Hashtable) htBBH.get(g1);
      for (int j=i+1; j<lGene.size(); j++) {
      	String g2 = lGene.get(j).toString();
      	Hashtable h2 = (Hashtable) htBBH.get(g2);
      	int nGenomes = 0;
      	double p = 0.0;
      	double p0 = 0.0;
     	
      	for (int k=0; k<lSpe.size(); k++) { /////////// start from 1 to not use query species
      	  if (h1.containsKey(lSpe.get(k)) && h2.containsKey(lSpe.get(k))) {
      	    Genome g = (Genome) htSpe.get(lSpe.get(k));
      	    String o1 = h1.get(lSpe.get(k)).toString();
      	    String o2 = h2.get(lSpe.get(k)).toString();
      	    
      	    if (g.containsGene(o1) && g.containsGene(o2)) {
	    if (k == 0) {
	      p0 = -Math.log(g.p[g.getOrderDist(o1, o2)]/0.5) / Math.log(10.0);
	    } else {
      	      nGenomes++;
      	      p += -Math.log(g.p[g.getOrderDist(o1, o2)]/0.5) / Math.log(10.0);
	    }
      	    }
      	  }	  
      	  
      	}
      	//System.out.println(g1+"\t"+g2+"\t"+p+"\t"+nGenomes+"\t"+p0);
      	double rounded = round24(p,8); 
       System.out.println(g1+"\t"+g2+"\t"+rounded);

      }      
    }
  }
  
  private static void readGenes(String fileGene) {
    Genome g = new Genome("", fileGene);
    htSpe.put("", g);
    lSpe.add("");
    
    try {
      BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileGene))));
      String sLine;
	while ((sLine = br.readLine()) != null) {
	  String[] sArray = sLine.split("\t");
	  lGene.add(sArray[0]);
	  Hashtable h = new Hashtable();
	  h.put("", sArray[0]);
	  htBBH.put(sArray[0], h);
	}
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
  
  private static void readSpe(String fileSpe) {
    try {
      BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileSpe))));
      String sLine;
	while ((sLine = br.readLine()) != null) {
	if (!sLine.equals(myspe)) {
	  htSpe.put(sLine, "");
	  lSpe.add(sLine);
	}
	}
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
  
  private static void readBBH(String fileBBH) {
    try {
      BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileBBH))));
      String sLine;
	while ((sLine = br.readLine()) != null) {
	  String[] sArray = sLine.split("\t");
	  String spe = sArray[1].substring(0,sArray[1].indexOf(':'));
	  if (htSpe.containsKey(spe)) {
	    //if (htBBH.containsKey(sArray[0])) {
	      Hashtable h = (Hashtable) htBBH.get(sArray[0]);
	      if (h != null)
	        h.put(spe, sArray[1]);
	    //} else {
	    //  Hashtable h = new Hashtable();
	    //  h.put(spe, sArray[1]);
	    //  htBBH.put(sArray[0], h);
	    //}
	  }
	}
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
  
  private static void readGenomes(String dir) {
    for (int i=1; i<lSpe.size(); i++) {
      String spe = lSpe.get(i).toString();
      String f = dir+spe+".out";
      Genome g = new Genome(spe, f);
      htSpe.put(spe, g);
    }
  }
  public static double round24(double d, int c) {
    int temp=(int)((d*Math.pow(10,c)));
    return (((double)temp)/Math.pow(10,c));
  }


}



class Genome {
  String name;
  Hashtable htGene2Ch;
  Hashtable htGeneOrder;
  int[] chSizes;
  int nGenes = 0;
  int maxSize = 0;
  double[] p;
  
  public Genome(String spe, String f) {
    name = spe;
    htGene2Ch = new Hashtable();
    htGeneOrder = new Hashtable();
    
    Hashtable htCh = new Hashtable();
    LinkedList lCh = new LinkedList();
    try {
      BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(f))));
      String sLine;
	while ((sLine = br.readLine()) != null) {
	  String[] sArray = sLine.split("\t");
	  htGene2Ch.put(sArray[0], sArray[1]);
	  Integer order = new Integer(sArray[6]);
	  htGeneOrder.put(sArray[0], order);
	  
	  if (!htCh.containsKey(sArray[1])) {
	    lCh.add(sArray[1]);
	    htCh.put(sArray[1], order);
	  }
	  if (order.compareTo((Integer) htCh.get(sArray[1]))>0)
	    htCh.put(sArray[1], order);
	}
    } catch (IOException e) {
      e.printStackTrace();
    }
    
    chSizes = new int[lCh.size()];
    for (int i=0; i<chSizes.length; i++) {
      Integer size = (Integer) htCh.get(lCh.get(i));
      chSizes[i] = size.intValue();
      nGenes += chSizes[i];
      maxSize = Math.max(chSizes[i], maxSize);
    }
    
    if (chSizes.length > 1) {	// linear
      p = new double[maxSize+2];
      for (int i=1; i<p.length; i++) {
      	for (int j=0; j<chSizes.length; j++) {
      	  if (i <= chSizes[j])
      	    p[i] += (i*(chSizes[j]-i)+i*(i-1.0)/2.0) / (nGenes*(nGenes-1)/2.0);
      	  else
      	    p[i] += chSizes[j]*1.0 / nGenes;
      	}
      }
    } else {	// circular
      p = new double[(maxSize+1)/2+1];
      for (int i=1; i<p.length; i++)
      	p[i] = Math.min(2.0*i/(maxSize-1), 1.0);
    }
    //System.out.println(name+"\t"+p.length+"\t"+p[p.length-1]);
  }
  
  public int getOrderDist(String g1, String g2) {
    if (htGene2Ch.get(g1).equals(htGene2Ch.get(g2))) {
      Integer o1 = (Integer) htGeneOrder.get(g1);
      Integer o2 = (Integer) htGeneOrder.get(g2);
      if (chSizes.length > 1) {	// linear
      	return Math.abs(o1.intValue()-o2.intValue());
      } else {	// circular
      	if (Math.abs(o1.intValue()-o2.intValue()) < 0.5*maxSize)
      	  return Math.abs(o1.intValue()-o2.intValue());
      	else
      	  return maxSize-Math.abs(o1.intValue()-o2.intValue());
      }
    }
    return maxSize+1; 
  }
  
  public boolean containsGene(String g) {
    return htGene2Ch.containsKey(g);
  }
  
}
