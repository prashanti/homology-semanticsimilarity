package experiment_3;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.io.RDFXMLOntologyFormat;
import org.semanticweb.owlapi.model.IRI;
import org.semanticweb.owlapi.model.OWLAxiom;
import org.semanticweb.owlapi.model.OWLClass;
import org.semanticweb.owlapi.model.OWLClassExpression;
import org.semanticweb.owlapi.model.OWLDataFactory;
import org.semanticweb.owlapi.model.OWLDeclarationAxiom;
import org.semanticweb.owlapi.model.OWLObjectProperty;
import org.semanticweb.owlapi.model.OWLOntology;
import org.semanticweb.owlapi.model.OWLOntologyCreationException;
import org.semanticweb.owlapi.model.OWLOntologyManager;
import org.semanticweb.owlapi.model.OWLOntologyStorageException;

public class AddWithHomologyGroupingClasses {
	public static void main(String[] args) throws OWLOntologyCreationException, OWLOntologyStorageException, IOException
	{
			
		String homologyontology="/Users/pmanda/Documents/phenotype-ontologies-read-only/src/ontology/" +
				"imports/uberon_import.owl";
		File file = new File(homologyontology);
		OWLOntologyManager manager = OWLManager.createOWLOntologyManager();
		OWLOntology ontology = manager.loadOntologyFromOntologyDocument(file);
		 String base ="http://purl.obolibrary.org/obo/";
		 
		 String ontologyname="/Users/pmanda/Documents/Github/HomologyAnalysis/data/Ontologies/" +
					"WithHomology.owl";
		 String outfile="/Users/pmanda/Documents/Github/" +
		    		"HomologyAnalysis/data/Experiment_3/WithHomologyGroupings.owl";
		File file1 = new File(ontologyname);
		OWLOntologyManager manager1 = OWLManager.createOWLOntologyManager();
		OWLOntology ontology1 = manager1.loadOntologyFromOntologyDocument(file1);
		OWLDataFactory factory = manager1.getOWLDataFactory();
		
		OWLObjectProperty inheresIn = factory.getOWLObjectProperty(IRI.create(base + "RO_0000052"));
		OWLObjectProperty homologousTo = factory.getOWLObjectProperty(IRI.create(base + "RO_HOM0000007"));
		OWLObjectProperty hasPart = factory.getOWLObjectProperty(IRI.create(base + "BFO_0000051"));
		OWLClass pato=factory.getOWLClass(IRI.create(base+ "PATO_0000001"));
		    
		String implicithomologyfile="/Users/pmanda/Documents/Github/HomologyAnalysis/data/homology-grouping-ids.txt";
		BufferedReader br = new BufferedReader(new FileReader(implicithomologyfile));
	    String line;
	    List<String> implicithomologyids = new ArrayList<String>();
	    
	   while ((line = br.readLine()) != null)
	   {
		   implicithomologyids.add(line.trim());
	   }
		   
		for (OWLClass uberonclass : ontology.getClassesInSignature(true))
		{
			String classname=""+uberonclass;
			if (!classname.contains("HOM") && !classname.contains("NCBI") 
			&& !classname.contains("Publication") )
			{
	    	if (classname.contains("http://purl.obolibrary.org/obo/"))
	    	{
	    	classname=classname.split("http://purl.obolibrary.org/obo/")[1];
	    	
	    	}
	    	classname=classname.replace(">","");
	    	classname=classname.replace("<","");
	    	System.out.println(classname);
	    	//if !Arrays.asList(implicithomologyids).contains(classname){
	    		
	    	// }

			// with homology class
			// has_part (Q and inheres_in homologous_to E)
	    	OWLClassExpression part1 = factory.getOWLObjectSomeValuesFrom(homologousTo,uberonclass); 
			OWLClassExpression part2 = factory.getOWLObjectSomeValuesFrom(inheresIn,part1);
			OWLClassExpression part3=factory.getOWLObjectIntersectionOf(pato,part2);
			OWLClassExpression withhomology = factory.getOWLObjectSomeValuesFrom(hasPart,part3);
			
			
			OWLClass NewClass1 = factory.getOWLClass(IRI.create(base + "WithHomology"+classname));
	        OWLDeclarationAxiom declarationAxiom = factory
	                .getOWLDeclarationAxiom(NewClass1);
	        manager.addAxiom(ontology1, declarationAxiom);
			OWLAxiom phenotypeEquivalence = factory.getOWLEquivalentClassesAxiom(NewClass1, withhomology);
	        manager1.addAxiom(ontology1, phenotypeEquivalence);
	        
	        // with homology class 
	        // has_part (Q and inheres_in (E or homologous_to E))
	        part1 = factory.getOWLObjectSomeValuesFrom(homologousTo,uberonclass); 
	    	part2 = factory.getOWLObjectUnionOf(uberonclass,part1);
			part3 = factory.getOWLObjectSomeValuesFrom(inheresIn,part2);
			OWLClassExpression part4=factory.getOWLObjectIntersectionOf(pato,part3);
			withhomology = factory.getOWLObjectSomeValuesFrom(hasPart,part4);
			
	        NewClass1 = factory.getOWLClass(IRI.create(base + "WithHomology2"+classname));
	        declarationAxiom = factory
	                .getOWLDeclarationAxiom(NewClass1);
	        manager.addAxiom(ontology1, declarationAxiom);
			phenotypeEquivalence = factory.getOWLEquivalentClassesAxiom(NewClass1, withhomology);
	        manager1.addAxiom(ontology1, phenotypeEquivalence);
	        
	        
	        
	        
	       // without homology class
	        // has_part (Q and inheres_in E)
	        part1 = factory.getOWLObjectSomeValuesFrom(inheresIn,uberonclass);
	        part2=factory.getOWLObjectIntersectionOf(pato,part1);
	        OWLClassExpression withouthomology = factory.getOWLObjectSomeValuesFrom(hasPart,part2);
	        OWLClass NewClass2 = factory.getOWLClass(IRI.create(base + "WithoutHomology"+classname));
	        OWLDeclarationAxiom declarationAxiom1 = factory
	                .getOWLDeclarationAxiom(NewClass2);
	        manager.addAxiom(ontology1, declarationAxiom1);
	        OWLAxiom phenotypeEquivalence2 = factory.getOWLEquivalentClassesAxiom(NewClass2, withouthomology);
	        manager1.addAxiom(ontology1, phenotypeEquivalence2);
	        
	        
	        		}
	}
		
		manager.saveOntology(ontology1,
                new RDFXMLOntologyFormat(),
                IRI.create((new File(outfile).toURI())));
    
	    	}
}
