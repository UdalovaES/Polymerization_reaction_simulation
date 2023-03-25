##Реакция полимеризации<a class="headerlink" href="#id1" title="Permalink to this heading"></a></h1>
<p>Простейшая реакция катализа записывается в виде: <span class="math notranslate nohighlight">\(A+B\xrightarrow{k}A+C\)</span>. Динамику концентрации вещества <span class="math notranslate nohighlight">\(A\)</span> и <span class="math notranslate nohighlight">\(B\)</span> при этом можно описать при помощи следующей системы дифференциальных уравнений:</p>
<blockquote>
<div><div class="math notranslate nohighlight">
\begin{split}\begin{cases}
\frac{dB}{dt} = -kAB,\\
\frac{dC}{dt} = kAB.
\end{cases}\end{split}</div>
</div></blockquote>
<p>В Системе уравнений, <span class="math notranslate nohighlight">\(A=\frac{N_{A}}{V_{0}}\)</span>, <span class="math notranslate nohighlight">\(B=\frac{N_{B}}{V_{0}}\)</span>, <span class="math notranslate nohighlight">\(C=\frac{N_{C}}{V_{0}}\)</span>, <span class="math notranslate nohighlight">\(k\)</span> — кинетическая константа скорости реакции, <span class="math notranslate nohighlight">\(N_{A}\)</span>, <span class="math notranslate nohighlight">\(N_{B}\)</span>, <span class="math notranslate nohighlight">\(N_{C}\)</span> — количество реагентов типа <span class="math notranslate nohighlight">\(A\)</span>, <span class="math notranslate nohighlight">\(B\)</span> и <span class="math notranslate nohighlight">\(C\)</span> соответственно, <span class="math notranslate nohighlight">\(V_{0}\)</span> — объем системы. Решение системы:</p>
<blockquote>
<div><div class="math notranslate nohighlight">
\begin{split}\begin{cases}
B = B_{0}e^{-kAt},\\
C = B_{0}-B=B_{0}(1-e^{-kAt}).
\end{cases}\end{split}</div>
</div></blockquote>
<p>Реагирующая система может быть описана в виде конечного объема со случайным начальным распределением <span class="math notranslate nohighlight">\(N\)</span> частиц в нем и заданными периодическими граничными условиями. Будем считать, что реакция происходит с некоторой вероятностью <span class="math notranslate nohighlight">\(p=0.001\)</span> когда частицы <span class="math notranslate nohighlight">\(A\)</span> и <span class="math notranslate nohighlight">\(B\)</span> находятся ближе чем <span class="math notranslate nohighlight">\(r_c=4.0\)</span> друг к другу. Взаимодействие частиц можно описать при помощи следущего потенциала:</p>
<blockquote>
<div><div class="math notranslate nohighlight">
\[V(\vec{r}_i,\vec{r}_j)=\varepsilon\left[A\left(\frac{\sigma}{r_{ij}}\right)^{4} + B\left(\frac{\sigma}{r_{ij}}\right)^2\right],\]</div>
</div></blockquote>
<p>где <span class="math notranslate nohighlight">\(\varepsilon=1\)</span> — энергетическая постоянная, <span class="math notranslate nohighlight">\(\sigma=1\)</span> — равновесное расстояние между частицами, <span class="math notranslate nohighlight">\(A\)</span> и <span class="math notranslate nohighlight">\(B\)</span> — константы, <span class="math notranslate nohighlight">\(r_{ij}=|\mathbf{r}_i-\mathbf{r}_j|\)</span> — текущее расстояние между частицами.  Если задать постоянные <span class="math notranslate nohighlight">\(A=1\)</span>, <span class="math notranslate nohighlight">\(B=-2\)</span>, то взаимодействие между частицами будет притягивающим, а при <span class="math notranslate nohighlight">\(A=0\)</span>, <span class="math notranslate nohighlight">\(B=1\)</span> — отталкивающим. Допустим, что продукт реакции <span class="math notranslate nohighlight">\(C\)</span> полимеризуется. Тогда результатом реакции будет полимер, состоящий из молекул типа <span class="math notranslate nohighlight">\(C\)</span>. Для описания процесса полимеризации предлагается использовать потенциал c притягивающими постоянными (<span class="math notranslate nohighlight">\(A=1\)</span> и <span class="math notranslate nohighlight">\(B=-2\)</span>) в случае, если обе частицы имеют тип <span class="math notranslate nohighlight">\(C\)</span> и отталкивающие (<span class="math notranslate nohighlight">\(A=0\)</span>, <span class="math notranslate nohighlight">\(B=1\)</span>) в другом случае.</p>
<p>Движение частиц описывается уравнениями Ланжевена:</p>
<blockquote>
<div><div class="math notranslate nohighlight">
\[m\frac{d^2\mathbf{r}_i}{dt^2}=-\nabla_iV-\lambda\frac{d\mathbf{r}_i}{dt}+\eta_i(t),\]</div>
</div></blockquote>
<p>где <span class="math notranslate nohighlight">\(m = 1.0\)</span> — масса частиц, <span class="math notranslate nohighlight">\(V({\mathbf{r}_i}) = V_{cov}({\mathbf{r}_i}) + V_{hic}({\mathbf{r}_i})\)</span> — потенциальная энергия системы, <span class="math notranslate nohighlight">\(\nabla_i V\)</span> — градиент потенциала по координате <span class="math notranslate nohighlight">\(i\)</span>-ой частицы. <span class="math notranslate nohighlight">\(\eta_i(t)\)</span> — случайная сила, описывающая соударения с молекулами воды с нормальным распределением:</p>
<blockquote>
<div><div class="math notranslate nohighlight">
\[\langle\eta_i(t)\eta_j(t)\rangle = 2\lambda k_BT\delta_{ij}\delta(t-t')\]</div>
</div></blockquote>
<p>Динамику Ланжевена можно рассматривать как динамику Ньютона с добавлением двух дополнительных сил: силы трения, пропорциональной скорости частицы и случайной силы, распределенной по Гауссу. Коэффициент пропорциональности <span class="math notranslate nohighlight">\(\lambda\)</span> называется коэфициентом демпфирования и задает силу влияния среды на молекулу. При больших <span class="math notranslate nohighlight">\(\lambda\)</span> система переходит к Брауновскому движению, при малых — к Ньютоновскому.</p>
<p>Для численного решения уравнения Ланжевена можно использовать следующую разностную схему:</p>
<blockquote>
<div><div class="math notranslate nohighlight">
\begin{align}\begin{aligned}\begin{split}\begin{split}
 m_i\frac{\mathbf{v}_{i}^{n+\frac{1}{2}}-\mathbf{v}_{i}^{n-\frac{1}{2}}}{\tau} &amp;= \mathbf{f}_{i}^{n} - \lambda m_i\frac{\mathbf{v}_{i}^{n+\frac{1}{2}}+\mathbf{v}_{i}^{n-\frac{1}{2}}}{2}+\sqrt{\frac{2k_BT\lambda m_i}{\tau}}\mathbf{r}_i^f\\
\frac{\mathbf{r}_{i}^{n+1}-\mathbf{r}_{i}^{n}}{\tau} &amp;= \mathbf{v}_{i}^{n+\frac{1}{2}}\end{split}\\\end{split}\end{aligned}\end{align}</div>
</div></blockquote>
<p>Здесь, <span class="math notranslate nohighlight">\(\lambda=50\)</span> — коэфициент демпфирования, <span class="math notranslate nohighlight">\(\mathbf{r}_i^f\)</span> — вектор из трех нормально распределенных случайных величин, с дисперсией <span class="math notranslate nohighlight">\(1\)</span> и математическим ожиданием <span class="math notranslate nohighlight">\(0\)</span>. Постоянная Больцмана <span class="math notranslate nohighlight">\(k_B=8.31\times10^{-3}\)</span> кДж/K*моль, шаг по времени <span class="math notranslate nohighlight">\(\tau=0.001\)</span> пс, температура <span class="math notranslate nohighlight">\(T=300\)</span> K.</p>
<section id="vmd">