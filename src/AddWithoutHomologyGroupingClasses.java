package experiment_3;

import java.io.File;
import java.io.FileNotFoundException;

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

public class AddWithoutHomologyGroupingClasses {
	public static void main(String[] args) throws OWLOntologyCreationException, FileNotFoundException, OWLOntologyStorageException
	{
			
		String homologyontology="/Users/pmanda/Documents/phenotype-ontologies-read-only/src/ontology/" +
				"imports/uberon_import.owl";
		File file = new File(homologyontology);
		OWLOntologyManager manager = OWLManager.createOWLOntologyManager();
		OWLOntology ontology = manager.loadOntologyFromOntologyDocument(file);
		 String base ="http://purl.obolibrary.org/obo/";
		 
		 String ontologyname="/Users/pmanda/Documents/Github/HomologyAnalysis/data/Ontologies/" +
					"WithoutHomology.owl";
		 String outfile="/Users/pmanda/Documents/Github/" +
		    		"HomologyAnalysis/data/Experiment_3/WithoutHomologyGroupings.owl";
		 
		File file1 = new File(ontologyname);
		OWLOntologyManager manager1 = OWLManager.createOWLOntologyManager();
		 
		 
			 OWLOntology ontology1 = manager1.loadOntologyFromOntologyDocument(file1);
				
		
		OWLDataFactory factory = manager1.getOWLDataFactory();
		 
		 
		 
		  OWLObjectProperty inheresIn = factory.getOWLObjectProperty(IRI.create(base + "RO_0000052"));
		    OWLObjectProperty hasPart = factory.getOWLObjectProperty(IRI.create(base + "BFO_0000051"));
		    OWLClass pato=factory.getOWLClass(IRI.create(base+ "PATO_0000001"));
		    
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
	    	// withouthomology class
	    	// has_part(Q and inheres_in E)
	    	
	        OWLClassExpression part1 = factory.getOWLObjectSomeValuesFrom(inheresIn,uberonclass);
	        OWLClassExpression part2=factory.getOWLObjectIntersectionOf(pato,part1);
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
