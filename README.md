# IACV-Registration-Project

**F12. Coarse point-clouds registration with 3D features**
			Given a pair of 3d point clouds in different reference frames, the aim of the project is to robustly
			estimate a roto-translation that coarsely registers them in a common reference frame. This can be done
			leveraging a RanSaC-like framework, where triplets of 3D features are used to estimate rorotranslations.
			To this end
			
			* a state-of-the-art 3D feature descriptor must be studied and implemented.
			*3D features must be matched
			*A least square method to instantiate roto-translations on triplets of matched features must be	implemented 
			*A RanSaC estimator must be implemented
