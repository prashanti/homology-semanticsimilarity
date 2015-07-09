package experiment_3;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Set;

import org.semanticweb.elk.owlapi.ElkReasonerFactory;
import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.model.OWLClass;
import org.semanticweb.owlapi.model.OWLOntology;
import org.semanticweb.owlapi.model.OWLOntologyCreationException;
import org.semanticweb.owlapi.model.OWLOntologyManager;
import org.semanticweb.owlapi.reasoner.NodeSet;
import org.semanticweb.owlapi.reasoner.OWLReasoner;
import org.semanticweb.owlapi.reasoner.OWLReasonerFactory;
public class getSuperClasses {



		public static void main(String[] args) throws OWLOntologyCreationException, FileNotFoundException
		{
			String ontologyname="/Users/pmanda/Documents/Github/" +
		    		"HomologyAnalysis/data/Experiment_3/" +
					"WithHomologyGroupings.owl";	
			File file = new File(ontologyname);
			OWLOntologyManager manager = OWLManager.createOWLOntologyManager();
			OWLOntology ontology = manager.loadOntologyFromOntologyDocument(file);
		    OWLReasonerFactory reasonerFactory = new ElkReasonerFactory();
		    OWLReasoner reasoner = reasonerFactory.createReasoner(ontology);
		    PrintWriter printWriter = new PrintWriter ("/Users/pmanda/Documents/Github/" +
		    		"HomologyAnalysis/data/Experiment_3/SubsumersWithHomology.txt");
		    for (OWLClass cls : ontology.getClassesInSignature(true))
		    {
		    	NodeSet<OWLClass> parents = reasoner.getSuperClasses(cls, false);
			    Set<OWLClass> clses = parents.getFlattened();
			    for (OWLClass parent : clses)
				    {
			    	 printWriter.println(cls+"\t"+parent);
				   }
			    
		    }
		    printWriter.close();
		    
		    
		    
			ontologyname="/Users/pmanda/Documents/Github/" +
		    		"HomologyAnalysis/data/Experiment_3/" +
					"WithoutHomologyGroupings.owl";	
			file = new File(ontologyname);
			 manager = OWLManager.createOWLOntologyManager();
			ontology = manager.loadOntologyFromOntologyDocument(file);
		     reasonerFactory = new ElkReasonerFactory();
		     reasoner = reasonerFactory.createReasoner(ontology);
		     printWriter = new PrintWriter ("/Users/pmanda/Documents/Github/" +
		    		"HomologyAnalysis/data/Experiment_3/SubsumersWithoutHomology.txt");
		    for (OWLClass cls : ontology.getClassesInSignature(true))
		    {
		    	NodeSet<OWLClass> parents = reasoner.getSuperClasses(cls, false);
			    Set<OWLClass> clses = parents.getFlattened();
			    for (OWLClass parent : clses)
				    {
			    	 printWriter.println(cls+"\t"+parent);
				    }
			    
		    }
		    printWriter.close();
	
		}
	}


