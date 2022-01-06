class ll_flt{

	private :

		typedef struct node{
			double info;
			node  *seg;
		} node;

		node *Topleft;
		node *newTopleft;
		node *Topright;

	public :

		// constructoras
		ll_flt();
		~ll_flt();
		void add(double info);
		void addrev(double info);


		// consultoras

		void resetpt(void **pt);
		void *noend(void *pt);
		double llpt(void *pt);
		void advpt(void **pt);

		// modificadoras

		void modpt(void *pt,double inf);

};
